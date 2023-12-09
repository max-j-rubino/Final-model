# WILD7400Finalproj
Wild 7400 final project code

This is a complex model following the ODD description

Class project: The purpose of this model is to evaluate the potential effects of the long-term harvest of the most vulnerable fish in the population. This model is inspired by the findings of Phillip et al. (2009). If angling removes the most vulnerable fish in the population, and if there is a genetic component of vulnerability we could expect high pressure angling to result in a population of fish with low angling vulnerability.

Key model attributes:

There is a heritable component to fish vulnerability. Fish with low vulnerability are less likely to be caught. Fish with high vulnerability are more likely to be caught and harvested.

There is an intercept in vulnerability based on a fish's habitat. Fish that are close to shore are more likely to be caught than fish that are far offshore. 

Fish grow following a Von Bertalanffy growth equation. Growth trajectory (K) is a heritable trait. 

Fishing and natural mortality happen at the same time.

Fish vulnerability to angling is bound between 0 and 1.

Purpose: The purpose of this model is to see if angling pressure could artificially select fish that are less vulnerable to angling. I hypothesize that under high angling pressure, the average vulnerability of fish will decrease as will catch per unit effort. Managing small impoundments often involves a high harvest of smaller fish. By removing smaller fish you are more likely to harvest fish with lower growth trajectories. Removal of these fish prevents them from contributing to the next generation of fish. Doing this should increase the mean length of fish at age. However, if vulnerability is heritable slow growing fish could evolve to have lower vulnerabilities. This would make the removal of small fish as a management strategy over time less effective. 

Overview:

Entities, state variables, and scales: This model has one entity which is fish. Each fish is described by three continuous variables: size (mm), growth trajectory, and vulnerability (0-1). Fish are also described by categorical scales, male or female, and nearshore or offshore.  The environment that the fish exist in is described as an additional two continuous variable natural mortality (M) and post-release mortality (a). Growth trajectories come from the Von Bertalanffy growth equations and the number of spawners per spawner. Spawners per spawner are lagged two years to account for time to mature. This model only uses adult fish. We assume juveniles are invulnerable to angling.

Processes: The probability of an individual surviving to spawn is the probability of a fish surviving to the time of spawning (S). S is defined as -exp(Z), where Z is the monthly survival rate. Z is equal to v[i](1-vr)+v[i]*vr*alpha+M. V[i] is the individual probability of capture. vr is the voluntary release rate, alpha is the post-release mortality rate, and M is the instantaneous monthly mortality rate. Length at age will be determined from a von Bertalanffy growth equation L(a)=Linf(1-exp(-k[i]((a-t0))) where Linf is the asymptotic length, k[i] is the growth coefficient for an individual fish, t0 is the theoretical age at which size is 0, and a is age. For 50 weeks fish grow and age weekly, and are acceptable to natural and angling mortality. After surviving 50 weeks each female will build a nest. Since males can fertilize more than one nest, they are paired with females randomly with replacement. Progeny vulnerability and growth coefficient are drawn from a beta distribution. The parameters of the beta distribution are drawn from the mean value of the parent's score for a trait and a predefined standard deviation.

Design Concepts:

Basic principles: This model describes concepts outlined in Phillip et al. (2009). If vulnerability to capture is heritable, we can expect angling to artificially select for fish less likely to be caught. This can lead to lower catch rates. We assume that there is a heritable component to each fish's probability of capture.

Emergence: The results of this model are the average vulnerability of fish in the population, the total mortality number of fish, the number of fish caught, the mean length and age of fish, and the number of trophy fish caught,  and catch per unit effort. Decreased catch per unit effort and decreased vulnerability demonstrate selective pressure towards less vulnerable individuals.

Adaptation: No direct objective seeking exists in this model. Indirect objective seeking exists that fish with higher vulnerability are less likely to reproduce. Fish that have low vulnerability are less likely to be caught and thus are more likely to have a higher contribution to the next generation.

Objectives: There is a tradeoff between the harvest selection of smaller fish, and the selection of lower vulnerability fish. At first, high harvest will likely be successful at removing the slower-growing individuals. But as the selective pressure begins to take effect, we could vulnerability decrease dramatically. There also may be a tradeoff between the harvest of smaller fish, and mortality of larger fish due to post-release morality. By having very high fishing pressure larger fish may experience substantial fishing mortality even when harvest rates are low. 

Learning: This model does not allow fish to make adaptive changes over time. A component could be added to allow anglers to make adaptive changes. Perhaps when catch rates are low they could leave.

Sensing: There is not a sensing component to this model.

Interactions: Fish interact when they reproduce, and when the carrying capacity of the pond begins to approach. As the population gets closer to K, reproductive success drops. 

Stochasticity: Stochasticity in this model results in the fish's weekly fates (caught and harvested, caught and released survivor, caught and released dead, not caught alive, not caught dead). These in turn result in whether or not the fish spawns or grows. Stochasitisty comes from these fates being drawn from a categorical distribution. Many of the probabilities in the categorical distribution are also derived from other distributions. For example, vulnerability is drawn from a beta distribution derived from the parents' mean vulnerability. There is also randomness in which individuals are paired to mate. Stochasticity is used to represent processes that are impossible to represent mechanically. These processes are either poorly understood or limited by computing power.

Collectives: Collectives are not represented in this model. However, there is density dependence in the number of eggs that survive to spawners, and the growth coefficient of fish.

Observations: Outputs from the model needed to observe population level change are: number of fish, population age and size structure, population size, number of fish caught, number of trophy fish caught, and mean vulnerability of fish.

Initialization: We assume equal sex ratios. Fish will be stocked at age 2 with a mean growth coefficient of 0.251 and sd of 0.002. Fishes vulnerability will be drawn from a BETA( distribution derived from a mean of 0.3 and sd of 0.05. 

Input Data: There is no input data for this model.

Subroutines: The model will have several subroutines. Weekly there will be the option for fish to be caught or not, and die or survive. They will also grow and age weekly. Fish cannot be caught twice in the same week. Fish that are not caught will have a greater survival, and be more likely to move to the next period. Yearly fish build nests and reproduce. The reproduction will result in progeny. These progeny will not enter the population until 2 years after reproduction to avoid having to model individual fish reproduction. All year 2 + fish can reproduce. Males are paired with females randomly with replacement. 
 
