function Groupe_28_partie2
% Groupe_28_partie2  organise le planning d'une usine de production de
% panneaux solaires. Cette fonction a ete realisee dans le cadre de la
% deuxieme partie du projet du cours de Modèles et méthodes d'optimisation
% I donne par Francois Glineur.
%
% Plus precisement, la fonction lit les donnees contenues dans les fichiers
% data_partie1.m et data_partie2.m. Ensuite, elle affiche les resultats
% demandes par les differentes questions de l'enonce et renvoie les
% resultats demandes par la Q13: gestion du personnel non relaxe.
%
% Q10: On minimise la meme fonction qu'en partie 1 et parmi les solutions
%      optimales, on trouve celles qui minimise et maximise les panneaux
%      sous-traites
%
% Q12: On introduit la gestion du personnel et on resoud le probleme relaxe
%
% Q13: On resoud le probleme non relaxe
%
% Entree: Le programme ne prend pas d'argument.
%
% Sorties: (de la Q13: gestion du personnel non relaxe)
%  plan_panneaux : vecteur T x 1. Contient le nombre de panneaux produits
%                  par l'usine chaque semaine.
%  plan_heure_sup : vecteur T x 1. Contient le nombre d'heures
%                   supplementaires totales faites chaque semaine.
%  plan_sous_traitant : vecteur T x 1. Contient le nombre de panneaux
%                       sous-traites chaque semaine.
%  plan_stockage : vecteur T x 1. Contient le nombre de panneaux mis en
%                  stock a la fin de chaque semaine. Ce vecteur correspond
%                  donc au nombre de panneaux stockes a la fin de la
%                  semaine, pas la production mis en stock et ajoutee au
%                  stock deja present.
%
%  plan_ouvriers : vecteur T x 1. Contient le nombre, entier, d'ouvriers
%                  que possede l'usine chaque semaine.
%
%  autre_vars : vecteur T x 5. Contient 5 vecteurs colonnes T x 1. Le
%               premier contient le nombre de panneaux fournis pour la
%               semaine i. Le deuxieme contient le nombre de panneaux
%               fournis pour la semaine i-1. Le troisieme contient le
%               nombre de panneaux fournis pour la semaine i-2. Le
%               quatrieme 
%
% Auteurs   : Andres Zarza Davila, Nizar Bousselmi et Gilles Peiffer
% Groupe    : 28

%Initialisation des donnees (peut eventuellement etre modifiees)
run('data_partie1.m');
run('data_partie2.m');
heures_par_jour=7;   %Heures de travail par ouvrier par jour
jours_par_semaine=5; %Jours de travail par ouvrier par semaine

%Variables (a ne pas modifier!)
heures_par_semaine=heures_par_jour*jours_par_semaine;
duree_assemblage=duree_assemblage/60;
salaire=T*nb_ouvriers*cout_horaire*nb_heures_remunerees;
I=ones(T,1);
E=speye(T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Q10: Probleme partie 1 avec seconde minimisation/maximisation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
fprintf('Q10: Minimisation des panneaux sous-traites\n\n')
fprintf('Nombre de panneaux sous-traites (minimise): %11.2f \n',sum(xMin(:,4)))

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
fprintf('\n\n')
  
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
fprintf('Q10: Maximisation des panneaux sous-traites\n\n')
fprintf('Nombre de panneaux sous-traites (maximise): %11.2f \n',sum(xMax(:,4)))

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
fprintf('\n\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Q12: Gestion du personnel relaxe %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
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

fprintf('Q12: Résolution du probleme relaxe\n\n')

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
fprintf('Production mise en stock  :%s\n\n',sprintf('%11.1f',mis_en_stock))
fprintf('Embauches/licencies       :%s\n',sprintf('%11.1f',x2(:,9)-x2(:,10)))
fprintf('Ouvriers                  :%s\n',sprintf('%11.1f',x2(:,8)))
fprintf('\n\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Q13: Gestion du personnel %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

optionsInt = optimoptions('intlinprog','Display','off');
intcon=(7*T+1):(8*T);

[xInt,fInt,exitFlagInt] = intlinprog(c2,intcon,A2,b2,Aeq2,beq2,lb2,ub2,optionsInt);
if exitFlagInt~=1
     disp('Le probleme ne semble pas avoir de solution (Q13) :-(')
     disp('La demande est probablement impossible a satisfaire avec les donnees fournies.')
     return
 end
xInt(intcon) = round(xInt(intcon));
xInt=reshape(xInt,[T,10]);

%Vecteurs a retourner
plan_panneaux=xInt(:,1);
plan_stockage=xInt(:,2);
plan_heure_sup=xInt(:,3);
plan_sous_traitant=xInt(:,4);
plan_stockage=xInt(:,2);
plan_ouvriers=xInt(:,8);
autres_vars=[xInt(:,5:7) xInt(:,9:10)];

%Affichage du tableau

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

fprintf('Q13: Résolution du probleme non-relaxe\n\n')

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
fprintf('Production mise en stock  :%s\n\n',sprintf('%11.1f',mis_en_stock))
fprintf('Embauches/licencies       :%s\n',sprintf('%11.1f',xInt(:,9)-xInt(:,10)))
fprintf('Ouvriers                  :%s\n',sprintf('%11.1f',xInt(:,8)))
fprintf('\n\n')

end