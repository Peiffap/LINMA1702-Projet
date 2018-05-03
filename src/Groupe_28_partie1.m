function x = Groupe_28_partie1(recalcul)
% Groupe_28_partie1  organise le planning d'une usine de production de
% panneaux solaires. Cette fonction a ete realisee dans le cadre du projet
% du cours de Modèles et méthodes d'optimisation I donne par Francois
% Glineur.
%
% Plus precisement, la fonction lit les donnees contenues dans le fichier
% data_partie1.m et renvoie les actions hebdomadaires que doit faire
% l'usine si elle veut minimiser son cout. Elle affiche dans un tableau les
% resultats de l'optimisation et dans un graphique une estimation de
% l'evolution du cout total si la demande est augmentee d'un vecteur
% epsilon*delta_demande ou epsilon est un scalaire entre 0 et 1. Ce sont
% les valeurs d'epsilon qui sont sur l'abscisse du graphique. Si le booleen
% "recalcul" est vrai (=1) alors le probleme est entierement recalcule avec
% les nouvelles valeurs de la demande et le nouveau cout exacte est affiche
% sur le graphique. Si recalcul=0 alors seulement l'estimation est
% calculee et affichee.
% Attention: si T est tres grand, le programme peut prendre un certain
% temps si on recalcule entierement le probleme pour chaque nouvelle
% demande.

% Entree:
%  recalcul : booleen. Si il est vrai (=1), en plus de l'estimation, le
%             probleme est entierement recalcule, et le nouveau cout total
%             affiche, pour chaque nouvelle valeur de
%             demande + epsilon*delta_demande. Si il est faux (=0),
%             seulement l'estimation est calculee et affichee.
%
% Sorties:
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
%  autre_vars : vecteur T x 3. Contient 3 vecteurs colonnes T x 1. Le
%               premier contient le nombre de panneaux fournis pour la
%               semaine i. Le deuxieme contient le nombre de panneaux
%               fournis pour la semaine i-1. Le troisieme contient le
%               nombre de panneaux fournis pour la semaine i-2.
%
% Auteurs   : Andres Zarza Davila, Nizar Bousselmi et Gilles Peiffer
% Groupe    : 28

%Initialisation des donnees (peut eventuellement etre modifiees)
run('data_partie1.m');
epsilon=0:.1:1;      %Vecteur epsilon de l'etude de la variation de la demande
heures_par_jour=7;   %Heures de travail par ouvrier par jour
jours_par_semaine=5; %Jours de travail par ouvrier par semaine

%Coming soon :-)
% cout_embauche=0;
% cout_licenciement=0;
% nb_max_ouvriers=0;

%Variables (a ne pas modifier!)
heures_par_semaine=heures_par_jour*jours_par_semaine;
duree_assemblage=duree_assemblage/60;
salaire=T*nb_ouvriers*cout_horaire*nb_heures_remunerees;
I=ones(T,1);
E=speye(T);
N=length(epsilon);
y=zeros(N,1);

%Choix de l'algorithme en fonction de la taille du probleme
if T>1000
    options = optimoptions('linprog','Display','off','Algorithm','interior-point');
else
    options = optimoptions('linprog','Display','off','Algorithm','dual-simplex');
end

%Verifie les donnees entrees
if recalcul ~=0 && recalcul~=1
    disp('Attention: recalcul doit valoir 1 si la valeur exacte de chaque nouvelle demande doit etre calculee,')
    disp('et 0 sinon et pas un autre nombre.')
end
if heures_par_jour>24 || heures_par_jour<0
    disp('ATTENTION: Nombre innatendu d''heures par jour :-(')
end
if jours_par_semaine>7 || jours_par_semaine<0
    disp('ATTENTION: Nombre innatendu de jours par semaine :-(')
end
if T<1
    disp('ATTENTION: Nombre innatendu de semaine :-(')
    disp('Le programme a ete arrete.')
    return
end
if nb_max_heure_sup<0
    fprintf('ATTENTION: nb_max_heure_sup= %d n''est pas consistant avec le fait que\n',nb_max_heure_sup)
    disp('les heures supplementaires doivent etre non-negatives.')
    disp('Le programme a ete arrete.')
    return
end
if nb_max_sous_traitant<0
    fprintf('ATTENTION: nb_max_sous_traitant= %d n''est pas consistant avec le fait que\n',nb_max_sous_traitant)
    disp('les panneaux achetes doivent etre non-negatifs.')
    disp('Le programme a ete arrete.')
    return
end

%Verifie que les ouvriers ne travaillent pas plus de 24h/jour
if heures_par_jour+nb_max_heure_sup >24
    disp('ATTENTION: il est possible que les ouvriers travaillent plus de 24h par jour avec ces donnees-ci.')
end

%Intialisation des matrices
c=[ I*cout_materiaux;
    I*cout_stockage;
    I*cout_heure_sup;
    I*cout_sous_traitant;
    sparse(T,1);
    I*penalite_1semaine;
    I*penalite_2semaines];

A=[ E*duree_assemblage sparse(T,T) -E sparse(T,4*T);
    sparse(2*T,2*T) speye(2*T) sparse(2*T,3*T)];

b=[ ones(T,1)*heures_par_semaine*nb_ouvriers;
    I*nb_max_heure_sup*nb_ouvriers;
    I*nb_max_sous_traitant];

if T==1
    Aeq=[ sparse(1,4*T) 1 sparse(1,2);
          1  -1 sparse(1,1) 1 -1 -1 -1];
    
    beq=[ demande';
          sparse(1,1)];
else
    Aeq=[ sparse(T,4*T) E sparse(T,1) E(:,1:end-1) sparse(T,2) E(:,1:end-2);
          E ([E(:,2:end) sparse(T,1)]-E) sparse(T,T) E -E -E -E;
          sparse(1,5*T) 1 sparse(1,2*T-1);  %x6,1=0
          sparse(1,6*T) 1 sparse(1,T-1);    %x7,1=0
          sparse(1,6*T+1) 1 sparse(1,T-2)]; %x7,2=0
    
    beq=[ demande';
          -stock_initial;
          sparse(T-2,1);
          stock_initial;
          sparse(3,1)]; %x6,1=x7,1=x7,2=0
end

lb=sparse(7*T,1);

%Cas initial
[x,objBase,exitFlag,~,dual] = linprog(c,A,b,Aeq,beq,lb,[],options);
if exitFlag~=1
    disp('Le probleme ne semble pas avoir de solution :-(')
    disp('La demande est probablement impossible a satisfaire avec les donnees fournies.')
    return
end

%Permet de passer d'un vecteur 7*T x 1 a une matrice T x 7
x = reshape(x,[T,7])

%Vecteurs a retourner
plan_panneaux=x(:,1);
plan_heure_sup=x(:,3);
plan_sous_traitant=x(:,4);
plan_stockage=x(:,2);
autres_vars=x(:,5:7);

objBase=objBase+salaire

%Stock initial de chaque semaine
stock_initial_semaine=[stock_initial plan_stockage(1:end-1)'];

%Production mise en stock
mis_en_stock=plan_stockage'-stock_initial_semaine;
mis_en_stock(end)=stock_initial-stock_initial_semaine(end);
tmp=mis_en_stock<0;
mis_en_stock(tmp)=0;

%Fournis apres 1 semaine
tmp1=x(:,6)';
fournis_1semaine=[tmp1(1,2:end) 0];

%Fournis apres 1 semaine
tmp3=x(:,7)';
fournis_2semaines=[tmp3(1,3:end) 0 0];
if T==1
    fournis_2semaines=0;
end

%Affichage du tableau
fprintf('Solution optimale calculee de valeur %11.2f \n',objBase)
fprintf('Semaines                  :%s\n',sprintf('%11d',1:T))
fprintf('Stock initial             :%s\n',sprintf('%11.1f',stock_initial_semaine))
fprintf('Production standard       :%s\n',sprintf('%11.1f',x(:,1)-x(:,3)/duree_assemblage))
fprintf('Production via heures sup :%s\n',sprintf('%11.1f',x(:,3)/duree_assemblage))
fprintf('Production sous-traitee   :%s\n',sprintf('%11.1f',x(:,4)))
fprintf('-> Production totale      :%s\n',sprintf('%11.1f',x(:,1)+x(:,4)))
fprintf('\n')
fprintf('Demande a satisfaire      :%s\n',sprintf('%11.1f',demande))
fprintf('Demande satis a temps     :%s\n',sprintf('%11.1f',x(:,5)))
fprintf('Demande non satisfaite    :%s\n',sprintf('%11.1f',demande'-x(:,5)))
fprintf('Demande satis apres 1 sem :%s\n',sprintf('%11.1f',fournis_1semaine))
fprintf('Demande satis apres 2 sem :%s\n',sprintf('%11.1f',fournis_2semaines))
fprintf('Production mise en stock  :%s\n',sprintf('%11.1f',mis_en_stock))
fprintf('\n')

%On recupere la solution du dual pour estimer la variation
solution_dual=dual.eqlin(1:T);
diff=-delta_demande*solution_dual;

%Parcourt des differentes valeurs de epsilon
%Ici on recalcule completement tout le probleme si recalcul=1
if recalcul==1
    for i=1:N
        %Modification du vecteur demande
        demande2=demande+epsilon(i)*delta_demande;
        if T==1
            beq=[ demande2';
                sparse(1,1)];
        else
            beq=[ demande2';
                -stock_initial;
                sparse(T-2,1);
                stock_initial;
                sparse(3,1)];
        end
        
        [~,obj,exitFlag] = linprog(c,A,b,Aeq,beq,lb,[],options);
        if exitFlag~=1
            fprintf('Le probleme ne semble pas avoir de solution en epsilon=%d (i=%d)\n',epsilon(i),i)
            disp('Le programme a ete arrete.')
            break
        else
            y(i)=obj;
        end
    end
    
    %Affichage du graphique
    figure
    plot(epsilon,y-objBase+salaire,'r','LineWidth',2)
    hold on
    p2=plot(epsilon,epsilon*diff,'--','LineWidth',2);
    p2.Color='b';
    title('Variation du cout total en fonction de la variation de la demande','FontSize',18)
    xlabel('epsilon','FontSize',16)
    ylabel('Variation du cout [euros]','FontSize',16)
    legend({'Valeur exacte','Estimation de la procedure'},'FontSize',15)
    
    %Si recalcul=0
else
    
    %Affichage du graphique
    figure
    p2=plot(epsilon,epsilon*diff,'--','LineWidth',2);
    p2.Color='b';
    title('Variation du cout total en fonction de la variation de la demande','FontSize',18)
    xlabel('epsilon','FontSize',16)
    ylabel('Variation du cout [euros]','FontSize',16)
    legend({'Estimation de la procedure'},'FontSize',15)
end
end