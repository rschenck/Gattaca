package GattacaExample;

import HAL.GridsAndAgents.AgentGrid0D;
import HAL.Rand;

import static HAL.Util.*;

import java.io.File;

class Parameters0D {

    public double mu=(1e-4);

    public double BIRTH_RATE = 0.1;
    public double DEATH_RATE = 0.099;
    public int n0 = 100;
    public int IGNORE_CLONES_BELOW_SIZE = 0;
    public boolean save_clonal_lineage = true;
    String[] AttributesList = new String[]{"Genome", "H", "S", "V"};
}

public class ExampleModel0d extends AgentGrid0D<Cell0D> {

    // tracking variables
    Parameters0D params = new Parameters0D();

    // neighborhoods
    int[]neighborhood=MooreHood(false);
    Rand rn;

    // initial gattaca parent genome
    Gattaca common_ancestor_genome;


    ExampleModel0d(int seed){
        super(Cell0D.class);

        rn = new Rand(seed);

        // does not appear anywhere in simulation: (0 popsize)
        common_ancestor_genome = new Gattaca(null, "", 0.5, 0.5, 0.3, rn);
        Gattaca clone1=new Gattaca(common_ancestor_genome, "",0.5, 0.5, 0.3, rn);

        int cells = 0;
        while (cells < params.n0) {

            Cell0D c = NewAgent().Init(clone1);
            c.genome.IncPop();

            cells++;
        }
    }

    // Step function for PD model ("steps" all cells through birth/death/mutation)
    void OriginalStep(){
        for (Cell0D c:this) {
            c.Step();
        }
        IncTick();
        CleanShuffle(rn);
    }


    /*
        Nonspatial()

        - run a simulation simulations: single non-spatial simulation
     */

    public static void Nonspatial(int totalTime, int modifier, boolean headless, int sims, boolean OVERWRITE) {

        String masterfoldername = "./GattacaExample/0d/";
        File dir = new File(masterfoldername);
        dir.mkdir();

        for (int seed = 0; seed < sims; seed++) {

            String foldername = masterfoldername;

            dir = new File(foldername);
            boolean success = dir.mkdir();

            if ((success) || OVERWRITE) {

                ExampleModel0d model = new ExampleModel0d(seed);
                Parameters0D p = model.params;

                for (int i = 0; i <= totalTime; i++) {

                    // step all models
                    model.OriginalStep();

                    // WRITE OUT VALUES
                    if (i % modifier == 0) {
                        System.out.println("Time: " + i + ", population: " + model.GetPop());


                        if (p.save_clonal_lineage) {
                            model.common_ancestor_genome.RecordClones(i);
                        }
                    }
                }

                // save clonal information in EvoFreq format:
                if (p.save_clonal_lineage) {
                    String filename = foldername + "gattaca_output"+ seed +".csv";
                    model.common_ancestor_genome.OutputClonesToCSV(filename, p.AttributesList, (Gattaca g) -> {
                        return GetAttributes(g);
                    }, p.IGNORE_CLONES_BELOW_SIZE);
                }
            }
        }

        return;
    }


    /*

        GetAttributes()
            - used for phylogeny tracker, to output clone-specific attributes

     */

    // Function to retrieve the attributes of your choice.
    public static String GetAttributes(Gattaca root) {
        return root.PrivateGenome + "," + Double.toString(root.h) + "," + Double.toString(root.s)+ "," + Double.toString(root.v);
    }

}