package GattacaExample;

import HAL.GridsAndAgents.AgentGrid2D;
import HAL.Gui.*;
import HAL.Rand;
import HAL.Util;

import static HAL.Util.*;

import java.io.File;

class Parameters2D {
    public double BIRTH_RATE = 0.1;
    public double DEATH_RATE = 0.05;
    public int sideLen = 500;
    public int n0 = sideLen*sideLen;
    public double IGNORE_CLONES_BELOW_FRACTION = 0.0; // (this is dependent on domain size)
    public boolean save_clonal_lineage = true;
    public int drawing_scaling_factor = 1;
    String[] AttributesList = new String[]{"Genome", "H", "S", "V"};
}

public class ExampleModel2d extends AgentGrid2D<Cell2D> {

    // tracking variables
    static final Parameters2D params = new Parameters2D();

    // neighborhoods
    int[]neighborhood=MooreHood(false);
    Rand rn;

    // initial gattaca parent genome
    Gattaca common_ancestor_genome;


    ExampleModel2d(int seed){
        super(params.sideLen,params.sideLen, Cell2D.class,false,false);

        rn = new Rand(seed);

        common_ancestor_genome = new Gattaca(null, "", 0.0, 0.0, 0.5, rn);
        Gattaca clone1=new Gattaca(common_ancestor_genome, "",0.0, 0.0, 0.5, rn);

        int cells = 0;
        while (cells < params.n0) {
            // random location
            int cell_id = rn.Int(this.length);
            Cell2D c  = this.GetAgent(cell_id);
            if (c==null) {
                NewAgentSQ(cell_id).Init(clone1);
                cells++;
            }

        }
    }

    // Step function for PD model ("steps" all cells through birth/death/mutation)
    void OriginalStep(){
        for (Cell2D c:this) {
            c.Step();
        }
        CleanShuffle(rn);
    }


    /*
        Square()

        - run a simulation simulations:
        - Single region with a homogeneous initial condition (e.g. ~50% density)
     */

    public static void Square(int totalTime, int modifier, boolean headless, int sims, int save_max, boolean OVERWRITE) {

        String masterfoldername = "./GattacaExample/2d/";
        File dir = new File(masterfoldername);
        dir.mkdir();

        for (int seed = 0; seed < sims; seed++) {

            String foldername = masterfoldername;

            dir = new File(foldername);
            boolean success = dir.mkdir();

            if ((success) || OVERWRITE) {


                ExampleModel2d model = new ExampleModel2d(seed);
                Parameters2D p = model.params;
                int delete_thresh = (int) (p.sideLen * p.sideLen * p.IGNORE_CLONES_BELOW_FRACTION);

                // VISUALIZE
                int vis_dim = p.sideLen + 4;
                UIGrid Vis = new UIGrid(vis_dim, vis_dim, p.drawing_scaling_factor);
                UIWindow win = (seed < save_max) ? CreateWindow(headless, Vis) : null;
                GifMaker smallGif = (seed < save_max) ? new GifMaker(foldername + "sim"+seed+".gif", 100, true) : null;

                for (int i = 0; i <= totalTime; i++) {

                    // step all models
                    model.OriginalStep();

                    // WRITE OUT VALUES
                    if (i % modifier == 0) {

                        System.out.println("Time: " + i + ", population: " + model.Pop());

                        Draw(model, Vis, i);

                        if (seed < save_max) {
                            smallGif.AddFrame(Vis);
                        }

                        if (p.save_clonal_lineage) {
                            model.common_ancestor_genome.RecordClones(i);
                        }
                    }
                }

                // save clonal information in EvoFreq format:
                if (p.save_clonal_lineage) {
                    if (seed < save_max) {
                        model.common_ancestor_genome.OutputClonesToCSV(foldername + "gattaca_output"+seed+".csv", p.AttributesList, (Gattaca g) -> {
                            return GetAttributes(g);
                        }, delete_thresh);
                    }
                }

                // CLOSE GIFS, WINDOWS
                if (seed < save_max) {
                    smallGif.Close();
                    if (!headless) {
                        win.Close();
                    }
                }
            }
        }

        return;
    }


    public static void Draw(ExampleModel2d model, UIGrid visCells, int time) {

        // color half by drivers and half by passengers
        for (int x = 0; x < model.params.sideLen; x++) {
            for (int y = 0; y < model.params.sideLen; y++) {

                Cell2D c = model.GetAgent(x,y);

                if (c== null) {
                    visCells.SetPix(x + 2,y + 2, Util.WHITE);
                } else {
                    visCells.SetPix(x + 2,y + 2,HSBColor(c.genome.h,c.genome.v,c.genome.s));
                }
            }
        }

//        visCells.SetString(Integer.toString((int) time), true, BLACK, WHITE, 2);

    }

    /*

        GetAttributes()
            - used for phylogeny tracker, to output clone-specific attributes

     */

    // Function to retrieve the attributes of your choice.
    public static String GetAttributes(Gattaca root) {
        return root.PrivateGenome + "," + Double.toString(root.h) + "," + Double.toString(root.s)+ "," + Double.toString(root.v);
    }


    public static UIWindow CreateWindow(boolean headless, UIGrid Vis) {
        UIWindow win = (headless) ? null : new UIWindow("Window", false, null, true);

        if (!headless) {
            win.AddCol(0, Vis);
            win.RunGui();
        }
        return win;
    }

}

