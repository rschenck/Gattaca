package GattacaExample;

import HAL.GridsAndAgents.Agent0D;

class Cell0D extends Agent0D<ExampleModel0d> {

    // each cell carries its own GattacaExample.Gattaca genome:
    Gattaca genome;

    Cell0D Init(Gattaca self_genome){
        genome = self_genome;
        if (G.params.save_clonal_lineage) {
            this.genome.IncPop();
        }
        return this;
    }
  
    Cell0D Divide(){

        Gattaca g = genome._RunPossibleMutation(G.GetTick());
        Cell0D c = G.NewAgent().Init(g);

        return c;
    }

    // constant death rate
    void Death() {
        if (G.params.save_clonal_lineage) {
            this.genome.DecPop();
        }
        Dispose();
    }

    void Step(){
        if(G.rn.Double()<(G.params.DEATH_RATE + G.params.BIRTH_RATE )){
            // death OR birth event occurs w/in this cell:
            double birth_event_probability = G.params.BIRTH_RATE/(G.params.DEATH_RATE + G.params.BIRTH_RATE);

            if(G.rn.Double()<(birth_event_probability )) {
                Divide();
            } else {
                Death();
                return;
            }
        }
    }
}

