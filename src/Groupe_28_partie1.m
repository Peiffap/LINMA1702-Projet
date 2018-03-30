%Initialisation des données
close all; clc;
run('data_partie1.m');
epsilon=0:.1:1;
heures_par_jour=7;
jours_par_semaine=5;

options = optimoptions('linprog','Display','off');
heures_par_semaine=heures_par_jour*jours_par_semaine;
duree_assemblage=duree_assemblage/60;
salaire=T*nb_ouvriers*cout_horaire*nb_heures_remunerees;
I=ones(T,1);
E=speye(T);
N=length(epsilon);
y=zeros(N,1);
    
%Vérifie T
if T<1
    disp('Nombre de semaine innatendu :-(')
    return
end

%Verifie que les ouvriers ne travaillent pas plus de 24h/d
if heures_par_jour+nb_max_heure_sup >24
    disp('Attention il est possible que les ouvriers travaillent plus de 24h par jour avec ces donnees-ci.')
end

%Indique si la solution sera entiere ou pas (pas au point)
if (35+nb_max_heure_sup)*nb_ouvriers/duree_assemblage==round((35+nb_max_heure_sup)*nb_ouvriers/duree_assemblage)
    disp('Entier')
else
    disp('Non-entier')
end

%Intialisation des matrices
disp('Version creuse')
tic

c=[ I*cout_materiaux;
    I*cout_stockage;
    I*cout_heure_sup;
    I*cout_sous_traitant;
    sparse(T,1);
    I*penalite_1semaine;
    I*penalite_2semaines];

A=[ E*duree_assemblage sparse(T,T) E*-nb_ouvriers sparse(T,4*T);
    sparse(2*T,2*T) speye(2*T) sparse(2*T,3*T)];

b=[ ones(T,1)*heures_par_semaine*nb_ouvriers;
    I*nb_max_heure_sup*nb_ouvriers;
    I*nb_max_sous_traitant];

if T==1
    Aeq=[sparse(1,4*T) 1 sparse(1,2);
        1 sparse(1,T-1) -1 sparse(1,2*T-1) 1 sparse(1,T-1) -1 sparse(1,T-1) -1 zeros(1,T-1) -1 zeros(1,T-1)];
    beq=[demande';
        sparse(1,1)];
else
    Aeq=[ sparse(T,4*T) E sparse(T,1) E(:,1:end-1) sparse(T,2) E(:,1:end-2);
        E ([E(:,2:end) sparse(T,1)]-E) sparse(T,T) E -E -E -E];
    
    beq=[ demande';
        -stock_initial;
        sparse(T-2,1);
        stock_initial];
end

lb=sparse(7*T,1);

%Cas initial
disp('Cas initial')
t=toc;
fprintf('Durée de l''initialisation des matrices: %d secondes\n',t)
tic;
[x,objBase,exitFlag,info,dual] = linprog(c,A,b,Aeq,beq,lb,[],options);
t=toc;
if exitFlag~=1
    disp('Le probleme ne semble pas avoir de solution.')
    return
end
fprintf('Durée de l''optimisation: %d secondes.\n',t)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% size(Aeq)
% AG=Aeq(1:T,:)
% j=1;
% for i=1:7*T
% 
%     if x(i)==0
%         AG(:,j)=[];
%     else
%     j=j+1;
%     end
% end
% size(AG)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = reshape(x,[T,7]);
objBase=objBase+salaire;
disp('Tableaux des resultats');
%disp([x demande'])
fprintf('Le cout total est %d euros.\n',objBase);
fprintf('La demande totale est de %d et nous fournissons %d dont %d produits et %d achetes.\n\n',sum(demande),sum(sum(x(:,5:7))),sum(sum(x(:,1))),sum(sum(x(:,4))));

solution_dual=dual.eqlin(1:T);
diff=-delta_demande*solution_dual;
disp('Comparaison des demandes')

%Parcourt des differentes valeurs de epsilon
for i=1:N
    
    demande=demande+epsilon(i)*delta_demande;
    
    if T==1
        beq=[demande';
            sparse(1,1)];
    else
        beq=[ demande';
            -stock_initial;
            sparse(T-2,1);
            stock_initial];
    end
    t=toc;
    fprintf('Iteration: %d (epsilon=%d)\n',i,epsilon(i))
    fprintf('Durée de l''initialisation des matrices: %d secondes\n',t)
    tic
    
    [x,obj,exitFlag] = linprog(c,A,b,Aeq,beq,lb,[],options);
    t=toc;
    if exitFlag~=1
        fprintf('Le probleme ne semble pas solvable en i=%d (epsilon=%d)\n',i,epsilon(i))
        break      
    else
        fprintf('Durée de l''optimisation: %d secondes.\n',t)
        y(i)=obj+salaire;
        fprintf('Augmentation: %d euros\n\n',y(i))
    end
end

%Affichage
figure
plot(epsilon,y-objBase,'r','LineWidth',2)
hold on
p2=plot(epsilon,epsilon*diff,'--','LineWidth',2);
p2.Color='b';
title('Variation du cout total en fonction de la variation de la demande')
xlabel('epsilon')
ylabel('Variation du cout [euros]')
legend('Valeur exacte','Estimation de la procedure')
