def Set_metapop(taillepop, nbsites): #Initialise une métapopulation de nbsites sites contenant chacun taillepop individus, un site aléatoire contiendra une individu infecté
    Metapop=[]

    for i in range(nbsites):
        if i == 0: Metapop.append([taillepop-1,1])
        else :Metapop.append([taillepop,0])
    Metapoparray = np.array(Metapop) # Les array sont stockés dans des blocs mémoire contigus et ca devrait nous accélérer entre 'un peu' et 'pas mal'
    return Metapoparray

def Get_xi(Metapop) : # Get density values of each subpoppulation of each site
    list_si = []
    list_ii = []
    for pop in Metapop :
        list_si.append(pop[0])
        list_ii.append(pop[1])
    list_xi = list_si+list_ii
    xi = np.array(list_xi)
    return xi

def GetMainMatrix(Metapop) : # Fonction qui récupère la matrice des propensities et celle qui contient les vecteurs de changement d'états

    # Séparation des S et I par convenance
    Slist = []
    Ilist = []

    for pop in Metapop :
        Slist.append(pop[0])
        Ilist.append(pop[1])
    Sarray = np.array(Slist)
    Iarray= np.array(Ilist)

    # Specifier la taille des matrices de sortie à l'avance
    Nbevents = 7
    NbPops = len(Sarray)
    # NbSpecies = len(Metapop) * len(Metapop[0])

    Nrow = Nbevents*NbPops
    Ncol = NbPops

    StateChangeMatrixS = np.zeros((Nrow, Ncol)) # Taille 70*10
    StateChangeMatrixI = np.zeros((Nrow, Ncol))
    MatrixPropensitiesS = np.zeros((Nrow, Ncol))
    MatrixPropensitiesI = np.zeros((Nrow, Ncol))

    # POur récupérer les ordres de réactions (utile après)
    Orders = []

    # print('Nombre de sites', NbPops)
    # print('Nombre entités', NbSpecies)
    # print('taille de matrice', StateChangeMatrixI.size)
    # print('Taille de boucle', range(len(Sarray)) )
    for specie in range(len(Sarray)) :
        propensities_I = np.zeros(NbPops)
        StatechangeVector_I = np.zeros(NbPops)

        propensities_S = np.zeros(NbPops)
        StatechangeVector_S = np.zeros(NbPops)

        # On rempli deux array de taille Nbspecie
        # L'un contient les propensity de chaque évènement
        # L'autre contient les changements d'états

        S = Sarray[specie]
        I = Iarray[specie]

        # Pour que les bons index soient remplis
        Index = Nbevents * specie # Va de Nbevent en Nbevent en partant de zéro

        # On décline tous les évènements possibles pour remplir les changements d'états associés + propensity
        # Reproduction S
        order = 1
        Orders.append(order)

        prop = r* S
        propensities_S[specie] = prop
        MatrixPropensitiesS[Index] = propensities_S

        changestate = +1
        StatechangeVector_S[specie] = changestate
        StateChangeMatrixS[Index] = StatechangeVector_S

        Index += 1

        # Mort S
        order = 3
        Orders.append(order)
        prop = r * S * (S+I) /k
        propensities_S[specie] = prop
        MatrixPropensitiesS[Index] = propensities_S

        changestate = -1
        StatechangeVector_S[specie] = changestate
        StateChangeMatrixS[Index] = StatechangeVector_S

        Index += 1

        # Migration S
        order = 1
        Orders.append(order)

        prop = d * S
        propensities_S[specie] = prop
        MatrixPropensitiesS[Index] = propensities_S

        changestate = -1
        StatechangeVector_S[specie] = changestate
        StateChangeMatrixS[Index] = StatechangeVector_S
        Index += 1
        # Migration I
        order = 1
        Orders.append(order)

        prop = d * I
        propensities_I[specie] = prop
        MatrixPropensitiesI[Index] = propensities_I

        changestate = -1
        StatechangeVector_I[specie] = changestate
        StateChangeMatrixI[Index]= StatechangeVector_I

        Index += 1
        # Mort I
        order = 1
        Orders.append(order)
        prop = alpha * I
        propensities_I[specie] = prop
        MatrixPropensitiesI[Index] = propensities_I

        changestate = -1
        StatechangeVector_I[specie] = changestate
        StateChangeMatrixI[Index]= StatechangeVector_I

        Index += 1
        # Guérison I
        order = 1
        Orders.append(order)

        prop = alpha * I
        propensities_I[specie] = prop
        MatrixPropensitiesI[Index] = propensities_I

        changestateS = +1
        changestateI = -1
        StatechangeVector_I[specie] = changestateI
        StateChangeMatrixI[Index]= StatechangeVector_I

        StatechangeVector_S[specie] = changestateS
        StateChangeMatrixS[Index] = StatechangeVector_S

        Index += 1
        # Infection
        order = 2
        Orders.append(order)

        prop = beta * S * I
        propensities_S[specie] = prop
        MatrixPropensitiesS[Index] = propensities_S

        changestateS = -1
        changestateI = +1

        StatechangeVector_I[specie] = changestateI
        StateChangeMatrixI[Index]= StatechangeVector_I

        StatechangeVector_S[specie] = changestateS
        StateChangeMatrixS[Index] = StatechangeVector_S

    StateChangeMatrix = np.concatenate((StateChangeMatrixS, StateChangeMatrixI), axis = 1)
    MatrixPropensities = np.concatenate((MatrixPropensitiesS, MatrixPropensitiesI), axis = 1)
    Reaction_orders = np.array(Orders)

    # Calcul de la somme des propensities
    aj = np.sum(MatrixPropensities, axis=1)
    return StateChangeMatrix, MatrixPropensities, Reaction_orders, aj

def GetCriticals(effectifs, propensities, statechange) : # Prends en entrée respectivement : Vecteur, Matrice, matrice
    ncrit = 11
    list_Indexcriticals = []
    dim_matrix = propensities.shape
    Critical_Matrix = np.zeros(dim_matrix)
    for index, effectif in enumerate(effectifs):
        if effectif < ncrit : list_Indexcriticals.append(index)
    for i in list_Indexcriticals : # Pour chaque pop critique

        for j in range(len(propensities)): # On regarde si aj > 0 ET vij < 0
            if propensities[j,i] > 0 and statechange[j,i] <0 : Critical_Matrix[j,i] =1
            else : pass
    vect_crit = np.sum(Critical_Matrix, axis= 1) # Axe 1 pour sommer sur les lignes !!!! Ne te gourres plus, ca va bien maintenant !
    return vect_crit # Renvoie un vecteur de booléens
