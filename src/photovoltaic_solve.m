function x = photovoltaic_solve()
    run('data_partie1.m')

    andres = 3;
    I = ones(T,1);
    E = eye(T);

    c=[ I * cout_materiaux;
        I * cout_stockage;
        I * nb_ouvriers * cout_heure_sup;
        I * cout_sous_traitant;
        zeros(T,1);
        I * penalite_1semaine;
        I * penalite_2semaines];

    %lb=zeros(7*T,1);
    lb = [zeros(2 * T - 1,1);
        stock_initial;
        zeros(5 * T,1)];
    ub = [I * Inf;
         ones(T - 1,1) * Inf;
         stock_initial;
         I * nb_max_heure_sup;
         I * nb_max_sous_traitant;
         I * Inf;
         I * Inf;
         I * Inf];

    A = [E * duree_assemblage / 60, zeros(T), E * -nb_ouvriers * 5,
        zeros(T,4 * T)];
    b = ones(T,1) * 35 * nb_ouvriers;

    Aeq = [E, -E, zeros(T), E, -E, -E, -E;
          zeros(T,4 * T), E, zeros(T,1), E(:,1:end - 1), zeros(T,2),
          E(:,1:end-2)];
    beq = [-stock_initial;
          zeros(T - 1,1);
          demande'];

    x = linprog(c, A, b, Aeq, beq, lb, ub);
    x = reshape(x, [T, 7]); %Vrai x

    %L = ["", "Produits", "Stockes", "Heures sup.", "Achetes", "Fournis pour i", "Fournis pour i-1", "Fournis pour i-2", "Demande"];
    %x = [(1:T)', x, demande']
    %x = [L; x];
end
