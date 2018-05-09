%Groupe_28_partie2
%Auteurs Nizar Andres Gilles
%Il faut ajouter les verifications de donnees
function Groupe_28_partie2
%Initialisation des donnees (peut eventuellement etre modifiees)
run('data_partie1.m');
heures_par_jour=7;   %Heures de travail par ouvrier par jour
jours_par_semaine=5; %Jours de travail par ouvrier par semaine

%Variables (a ne pas modifier!)
heures_par_semaine=heures_par_jour*jours_par_semaine;
duree_assemblage=duree_assemblage/60;
salaire=T*nb_ouvriers*cout_horaire*nb_heures_remunerees;
I=ones(T,1);
E=speye(T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Q10:Probleme partie 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Choix de l'algorithme en fonction de la taille du probleme
%n+m = 12*T+3 > 10 000 => interior-point
if T>1000
    options = optimoptions('linprog','Display','off','Algorithm','interior-point');
else
    options = optimoptions('linprog','Display','off','Algorithm','dual-simplex');
end

%Intialisation des matrices
c1=[ I*cout_materiaux;
     I*cout_stockage;
     I*cout_heure_sup;
     I*cout_sous_traitant;
     sparse(T,1);
     I*penalite_1semaine;
     I*penalite_2semaines];

A=[ E*duree_assemblage sparse(T,T) -E sparse(T,4*T);
    sparse(T,2*T) E sparse(T,4*T)];

b=[ ones(T,1)*heures_par_semaine*nb_ouvriers;
    I*nb_max_heure_sup*nb_ouvriers];

if T==1
    Aeq=[ sparse(1,4) 1 sparse(1,2);
          1 -1 sparse(1,1) 1 -1 -1 -1];
    
    beq=[ demande';
          sparse(1,1)];
    
    ub=[ inf;
         inf;
         inf;
         nb_max_sous_traitant;
         inf;
         sparse(1,1);  %x6,1=0
         sparse(1,1)]; %x7,1=0
    
else
    
    Aeq=[ sparse(T,4*T) E sparse(T,1) E(:,1:end-1) sparse(T,2) E(:,1:end-2);
          E ([E(:,2:end) sparse(T,1)]-E) sparse(T,T) E -E -E -E];
    
    beq=[ demande';
          -stock_initial;
          sparse(T-2,1);
          stock_initial];
    
    ub=[ ones(3*T,1)*inf;
         I*nb_max_sous_traitant;
         I*inf;
         sparse(1,1);  %x6,1=0
         ones(T-1,1)*inf;
         sparse(2,1);  %x7,1=x7,2=0
         ones(T-2,1)*inf];
end

lb=sparse(7*T,1);

[~,f,exitFlag] = linprog(c1,A,b,Aeq,beq,lb,ub,options);
if exitFlag~=1
    disp('Le probleme ne semble pas avoir de solution (Q10) :-(')
    disp('La demande est probablement impossible a satisfaire avec les donnees fournies.')
    return
end

%Nouvelle contrainte
AdoubleOpt=[Aeq;
            c1'];
bdoubleOpt=[beq;
            f];
        
%Minimisation des panneaux sous-traites
cMin=[ sparse(3*T,1);
       I*cout_sous_traitant;
       sparse(3*T,1)];
xMin = linprog(cMin,A,b,AdoubleOpt,bdoubleOpt,lb,ub,options);
xMin = reshape(xMin,[T,7]);
 
%Affichage du tableau

plan_stockage=xMin(:,2);
%Stock initial de chaque semaine
stock_initial_semaine=[stock_initial plan_stockage(1:end-1)'];

%Production mise en stock
mis_en_stock=plan_stockage'-stock_initial_semaine;
mis_en_stock(end)=stock_initial-stock_initial_semaine(end);
tmp=mis_en_stock<0;
mis_en_stock(tmp)=0;

%Fournis apres 1 semaine
tmp1=xMin(:,6)';
fournis_1semaine=[tmp1(1,2:end) 0];

%Fournis apres 1 semaine
tmp3=xMin(:,7)';
fournis_2semaines=[tmp3(1,3:end) 0 0];
if T==1
    fournis_2semaines=0;
end
fprintf('Q10: Minimisation des panneaux sous-traites\n')
fprintf('Cout des panneaux sous-traites (minimise): %11.2f \n',sum(xMin(:,4)))

fprintf('Solution optimale calculee de valeur %11.2f \n',f+salaire)
fprintf('Semaines                  :%s\n',sprintf('%11d',1:T))
fprintf('Stock initial             :%s\n',sprintf('%11.1f',stock_initial_semaine))
fprintf('Production standard       :%s\n',sprintf('%11.1f',xMin(:,1)-xMin(:,3)/duree_assemblage))
fprintf('Production via heures sup :%s\n',sprintf('%11.1f',xMin(:,3)/duree_assemblage))
fprintf('Production sous-traitee   :%s\n',sprintf('%11.1f',xMin(:,4)))
fprintf('-> Production totale      :%s\n',sprintf('%11.1f',xMin(:,1)+xMin(:,4)))
fprintf('\n')
fprintf('Demande a satisfaire      :%s\n',sprintf('%11.1f',demande))
fprintf('Demande satis a temps     :%s\n',sprintf('%11.1f',xMin(:,5)))
fprintf('Demande non satisfaite    :%s\n',sprintf('%11.1f',demande'-xMin(:,5)))
fprintf('Demande satis apres 1 sem :%s\n',sprintf('%11.1f',fournis_1semaine))
fprintf('Demande satis apres 2 sem :%s\n',sprintf('%11.1f',fournis_2semaines))
fprintf('Production mise en stock  :%s\n',sprintf('%11.1f',mis_en_stock))
fprintf('\n')
  
%Maximisation des panneaux sous-traites
xMax = linprog(-cMin,A,b,AdoubleOpt,bdoubleOpt,lb,ub,options);
xMax = reshape(xMax,[T,7]);

%Affichage du tableau

plan_stockage=xMax(:,2);
%Stock initial de chaque semaine
stock_initial_semaine=[stock_initial plan_stockage(1:end-1)'];

%Production mise en stock
mis_en_stock=plan_stockage'-stock_initial_semaine;
mis_en_stock(end)=stock_initial-stock_initial_semaine(end);
tmp=mis_en_stock<0;
mis_en_stock(tmp)=0;

%Fournis apres 1 semaine
tmp1=xMax(:,6)';
fournis_1semaine=[tmp1(1,2:end) 0];

%Fournis apres 1 semaine
tmp3=xMax(:,7)';
fournis_2semaines=[tmp3(1,3:end) 0 0];
if T==1
    fournis_2semaines=0;
end
fprintf('Q10: Maximisation des panneaux sous-traites\n')
fprintf('Cout des panneaux sous-traites (maximise): %11.2f \n',sum(xMax(:,4)))

fprintf('Solution optimale calculee de valeur %11.2f \n',f+salaire)
fprintf('Semaines                  :%s\n',sprintf('%11d',1:T))
fprintf('Stock initial             :%s\n',sprintf('%11.1f',stock_initial_semaine))
fprintf('Production standard       :%s\n',sprintf('%11.1f',xMax(:,1)-xMax(:,3)/duree_assemblage))
fprintf('Production via heures sup :%s\n',sprintf('%11.1f',xMax(:,3)/duree_assemblage))
fprintf('Production sous-traitee   :%s\n',sprintf('%11.1f',xMax(:,4)))
fprintf('-> Production totale      :%s\n',sprintf('%11.1f',xMax(:,1)+xMax(:,4)))
fprintf('\n')
fprintf('Demande a satisfaire      :%s\n',sprintf('%11.1f',demande))
fprintf('Demande satis a temps     :%s\n',sprintf('%11.1f',xMax(:,5)))
fprintf('Demande non satisfaite    :%s\n',sprintf('%11.1f',demande'-xMax(:,5)))
fprintf('Demande satis apres 1 sem :%s\n',sprintf('%11.1f',fournis_1semaine))
fprintf('Demande satis apres 2 sem :%s\n',sprintf('%11.1f',fournis_2semaines))
fprintf('Production mise en stock  :%s\n',sprintf('%11.1f',mis_en_stock))
fprintf('\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Q12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
c2=[ c1;
     I*cout_horaire*nb_heures_remunerees;
     I*cout_embauche;
     I*cout_licenciement];
 
A2=[ E*duree_assemblage sparse(T,T) -E sparse(T,4*T) -heures_par_semaine*E sparse(T,2*T);
     sparse(T,2*T) E sparse(T,4*T) -E*nb_max_heure_sup sparse(T,2*T)];

b2=sparse(2*T,1);

if T==1
    Aeq2=[ sparse(1,4) 1 sparse(1,5);
           1 -1 sparse(1,1) 1 -1 -1 -1 sparse(1,3);
           sparse(1,7) 1 -1 1];
    
    beq2=[ demande';
           sparse(1,1);
           nb_ouvriers];
    
    ub2=[ inf;
          inf;
          inf;
          nb_max_sous_traitant;
          inf;
          sparse(1,1);  %x6,1=0
          sparse(1,1);  %x7,1=0
          nb_max_ouvriers;
          inf;
          inf];
else
    Aeq2=[ [Aeq sparse(2*T,3*T)];
           sparse(T,7*T) (E-[E(:,2:end) sparse(T,1)]) -E E];
    
    beq2=[ beq;
           nb_ouvriers;
           sparse(T-1,1)];
    
    ub2=[ ub;
          I*nb_max_ouvriers;
          ones(2*T,1)*inf];
    
end

lb2=sparse(10*T,1);

%Choix de l'algorithme en fonction de la taille du probleme
%n+m = 17*T+3 > 10 000 => interior-point
if T>600
    options2 = optimoptions('linprog','Display','off','Algorithm','interior-point');
else
    options2 = optimoptions('linprog','Display','off','Algorithm','dual-simplex');
end
    
[x2,f2,exitFlag] = linprog(c2,A2,b2,Aeq2,beq2,lb2,ub2,options2);
if exitFlag~=1
    disp('Le probleme ne semble pas avoir de solution (Q12) :-(')
    disp('La demande est probablement impossible a satisfaire avec les donnees fournies.')
    return
end
x2=reshape(x2,[T,10]);

%Affichage du tableau

plan_stockage=x2(:,2);
%Stock initial de chaque semaine
stock_initial_semaine=[stock_initial plan_stockage(1:end-1)'];

%Production mise en stock
mis_en_stock=plan_stockage'-stock_initial_semaine;
mis_en_stock(end)=stock_initial-stock_initial_semaine(end);
tmp=mis_en_stock<0;
mis_en_stock(tmp)=0;

%Fournis apres 1 semaine
tmp1=x2(:,6)';
fournis_1semaine=[tmp1(1,2:end) 0];

%Fournis apres 1 semaine
tmp3=x2(:,7)';
fournis_2semaines=[tmp3(1,3:end) 0 0];
if T==1
    fournis_2semaines=0;
end

fprintf('Q12: Résolution du probleme relaxe\n')
fprintf('Cout des panneaux sous-traites (maximise): %11.2f \n',sum(x2(:,4)))

fprintf('Solution optimale calculee de valeur %11.2f \n',f2)
fprintf('Semaines                  :%s\n',sprintf('%11d',1:T))
fprintf('Stock initial             :%s\n',sprintf('%11.1f',stock_initial_semaine))
fprintf('Production standard       :%s\n',sprintf('%11.1f',x2(:,1)-x2(:,3)/duree_assemblage))
fprintf('Production via heures sup :%s\n',sprintf('%11.1f',x2(:,3)/duree_assemblage))
fprintf('Production sous-traitee   :%s\n',sprintf('%11.1f',x2(:,4)))
fprintf('-> Production totale      :%s\n',sprintf('%11.1f',x2(:,1)+x2(:,4)))
fprintf('\n')
fprintf('Demande a satisfaire      :%s\n',sprintf('%11.1f',demande))
fprintf('Demande satis a temps     :%s\n',sprintf('%11.1f',x2(:,5)))
fprintf('Demande non satisfaite    :%s\n',sprintf('%11.1f',demande'-x2(:,5)))
fprintf('Demande satis apres 1 sem :%s\n',sprintf('%11.1f',fournis_1semaine))
fprintf('Demande satis apres 2 sem :%s\n',sprintf('%11.1f',fournis_2semaines))
fprintf('Production mise en stock  :%s\n',sprintf('%11.1f',mis_en_stock))
fprintf('Embauches/licencies       :%s\n',sprintf('%11.1f',x2(:,9)-x2(:,10)))
fprintf('Ouvriers                  :%s\n',sprintf('%11.1f',x2(:,8)))
fprintf('\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Q13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

optionsInt = optimoptions('intlinprog','Display','off');
intcon=(7*T+1):(10*T);

[xInt,fInt,exitFlagInt] = intlinprog(c2,intcon,A2,b2,Aeq2,beq2,lb2,ub2,optionsInt);
if exitFlagInt~=1
     disp('Le probleme ne semble pas avoir de solution (Q13) :-(')
     disp('La demande est probablement impossible a satisfaire avec les donnees fournies.')
     return
 end
xInt(intcon) = round(xInt(intcon));
xInt=reshape(xInt,[T,10]);
%Affichage du tableau

plan_stockage=xInt(:,2);
%Stock initial de chaque semaine
stock_initial_semaine=[stock_initial plan_stockage(1:end-1)'];

%Production mise en stock
mis_en_stock=plan_stockage'-stock_initial_semaine;
mis_en_stock(end)=stock_initial-stock_initial_semaine(end);
tmp=mis_en_stock<0;
mis_en_stock(tmp)=0;

%Fournis apres 1 semaine
tmp1=xInt(:,6)';
fournis_1semaine=[tmp1(1,2:end) 0];

%Fournis apres 1 semaine
tmp3=xInt(:,7)';
fournis_2semaines=[tmp3(1,3:end) 0 0];
if T==1
    fournis_2semaines=0;
end

fprintf('Q13: Résolution du probleme non-relaxe\n')
fprintf('Cout des panneaux sous-traites (maximise): %11.2f \n',sum(xInt(:,4)))

fprintf('Solution optimale calculee de valeur %11.2f \n',fInt)
fprintf('Semaines                  :%s\n',sprintf('%11d',1:T))
fprintf('Stock initial             :%s\n',sprintf('%11.1f',stock_initial_semaine))
fprintf('Production standard       :%s\n',sprintf('%11.1f',xInt(:,1)-xInt(:,3)/duree_assemblage))
fprintf('Production via heures sup :%s\n',sprintf('%11.1f',xInt(:,3)/duree_assemblage))
fprintf('Production sous-traitee   :%s\n',sprintf('%11.1f',xInt(:,4)))
fprintf('-> Production totale      :%s\n',sprintf('%11.1f',xInt(:,1)+x2(:,4)))
fprintf('\n')
fprintf('Demande a satisfaire      :%s\n',sprintf('%11.1f',demande))
fprintf('Demande satis a temps     :%s\n',sprintf('%11.1f',xInt(:,5)))
fprintf('Demande non satisfaite    :%s\n',sprintf('%11.1f',demande'-xInt(:,5)))
fprintf('Demande satis apres 1 sem :%s\n',sprintf('%11.1f',fournis_1semaine))
fprintf('Demande satis apres 2 sem :%s\n',sprintf('%11.1f',fournis_2semaines))
fprintf('Production mise en stock  :%s\n',sprintf('%11.1f',mis_en_stock))
fprintf('Embauches/licencies       :%s\n',sprintf('%11.1f',xInt(:,9)-xInt(:,10)))
fprintf('Ouvriers                  :%s\n',sprintf('%11.1f',xInt(:,8)))
fprintf('\n')

%Inutile pour l'instant
% 
% %Permet de passer d'un vecteur 7*T x 1 a une matrice T x 7
% x = reshape(x,[T,10]);
% 
% %Vecteurs a retourner
% plan_panneaux=x(:,1);
% plan_heure_sup=x(:,3);
% plan_sous_traitant=x(:,4);
% plan_stockage=x(:,2);
% autres_vars=x(:,5:7);

end