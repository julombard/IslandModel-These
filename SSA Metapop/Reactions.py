
############# Définition des réactions #####################

# Reproduction S
#Rate et Stoechiometrie
r*S
S += 1

# Mort S
r*S*(S+I)/k
S -= 1

#Infection
beta * S * I
I += 1
S-= 1

Migration S
d * S
S-=1
#Sélection du site receveur
Site = random.randint(1, len(Metapop))
Metapop[Site][0] += 1

#Guerison:
gamma * I
I+= 1
S -= 1

#Mort I
alpha * I
I -= 1

# Migration I
d* I
I -= I
#Sélection du site receveur
Site = random.randint(1, len(Metapop))
Metapop[Site][1] += 1




