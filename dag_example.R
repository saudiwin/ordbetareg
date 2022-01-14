library(dagitty)
library(ggdag)

# roommate causal graph

roommate <- dagitty("dag{
        d -> r;
        p -> r;
        s -> r;
        n -> r;
        t -> r;
        v -> r;
        l -> r;
        d -> l;
        p  -> l;
        p ->  n;
        v -> s;
        v -> n;
        l -> v;
}")

ggdag(roommate)

# find variables that aren't related
impliedConditionalIndependencies(roommate)

# find the adjustment set
adjustmentSets(roommate, exposure="l", outcome =  "r")

# ggdag version
ggdag_adjustment_set(roommate, exposure="l",outcome="r")
