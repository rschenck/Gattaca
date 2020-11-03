package GattacaExample;

import HAL.GridsAndAgents.AgentSQ3Dunstackable;

class Cell3D extends AgentSQ3Dunstackable<ExampleModel3d> {

    // each cell carries its own GattacaExample.Gattaca genome:
    Gattaca genome;

    Cell3D Init(Gattaca self_genome){
        genome = self_genome;
        if (G.params.save_clonal_lineage) {
            this.genome.IncPop();
        }
        return this;
    }

    Cell3D Mutate(){
        if (G.rn.Double() <  G.params.mu) {
            // initiate new clone, with random color:
            this.genome.DecPop();
            this.genome = new Gattaca(this.genome, "", G.rn.Double(), G.rn.Double(), G.rn.Double(), G.rn);
            this.genome.IncPop();
        }
        return this;
    }

    Cell3D Divide(){
        int nDivOptions = G.MapEmptyHood(G.neighborhood,Xsq(),Ysq(),Zsq());
        if(nDivOptions==0){ return null; }
        int nextAgentID = G.neighborhood[G.rn.Int(nDivOptions)];
        return G.NewAgentSQ(nextAgentID).Init(genome).Mutate();
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