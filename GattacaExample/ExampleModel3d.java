package GattacaExample;
import HAL.GridsAndAgents.AgentGrid3D;
import HAL.Gui.*;
import HAL.Rand;
import static HAL.Util.*;
import HAL.Gui.OpenGL3DWindow;
import HAL.Util;

import java.io.File;

class Parameters3D {
    public double mu=(1e-4);
    public double BIRTH_RATE = 0.1;
    public double DEATH_RATE = 0.099;
    public int sideLen = 50;
    public int n0 = sideLen*sideLen;
    public double IGNORE_CLONES_BELOW_FRACTION = 0.0; // (this is dependent on domain size)
    public boolean save_clonal_lineage = true;
    public int drawing_scaling_factor = 3;
    String[] AttributesList = new String[]{"Genome", "H", "S", "V"};

}

public class ExampleModel3d extends AgentGrid3D<Cell3D> {

    // tracking variables
    static final Parameters3D params = new Parameters3D();

    // initial gattaca parent genome
    Gattaca common_ancestor_genome;

    int[]neighborhood = MooreHood3D(false);
    Rand rn;

    ExampleModel3d(int seed){
        super(params.sideLen,params.sideLen, params.sideLen, Cell3D.class, false, false, false);

        rn = new Rand(seed);

        common_ancestor_genome = new Gattaca(null, "", 0.5, 0.5, 0.3, rn);

        int cells = 0;
        while (cells < params.n0) {
            // random location
            int cell_id = rn.Int(this.length);
            Cell3D c  = this.GetAgent(cell_id);
            if (c==null) {
                NewAgentSQ(cell_id).Init(common_ancestor_genome);
                cells++;
            }
        }

    }

    void OriginalStep(){
        for (Cell3D c:this) {
            c.Step();
        }
        CleanShuffle(rn);
    }


    /*
        Cube()

            - runs a single simulation for parameters3D.totalTime time steps
            - invasion occurs after <kd> = 2
                - this happens by setting bigModel.p.INVASION = true;
            - outputs Muller information, basic population statistics ("everything" arrays) and GIFs
     */

    public static void Cube(int totalTime, int modifier, boolean headless, int sims, int save_max, boolean OVERWRITE) {


        String masterfoldername = "3d/";
        File dir = new File(masterfoldername);
        dir.mkdir();

        for (int seed = 0; seed < sims; seed++) {

            String foldername = masterfoldername; // + "seed" + Integer.toString((int) (seed)) + "/";

            dir = new File(foldername);
            boolean success = dir.mkdir();

            if ((success) || OVERWRITE) {


                ExampleModel3d model = new ExampleModel3d(seed);

                Parameters3D p = model.params;

                int delete_thresh = (int) (p.sideLen * p.sideLen * p.IGNORE_CLONES_BELOW_FRACTION);

                // VISUALIZE
                int vis_dim = p.sideLen + 4;

                System.out.println(vis_dim);

                UIGrid Vis = new UIGrid(vis_dim, vis_dim, p.drawing_scaling_factor);
                UIWindow win = (seed < save_max) ? CreateWindow(headless, Vis) : null;
                GifMaker gif = (seed < save_max) ? new GifMaker(foldername + "sim.gif", 100, true) : null;



                for (int i = 0; i <= totalTime; i++) {

                    // step all models
                    model.OriginalStep();

                    // WRITE OUT VALUES
                    if (i % modifier == 0) {

                        DrawProjection(Vis, p, model,i);

                        if (seed < save_max) {
                            gif.AddFrame(Vis);
                        }

                        if (p.save_clonal_lineage) {
                            model.common_ancestor_genome.RecordClones(i);
                        }
                    }
                }

                // save clonal information in EvoFreq format:
                if (p.save_clonal_lineage) {
                    if (seed < save_max) {
                        model.common_ancestor_genome.OutputClonesToCSV(foldername + "gattaca_output.csv", p.AttributesList, (Gattaca g) -> {
                            return GetAttributes(g);
                        }, delete_thresh);
                    }
                }

                // CLOSE GIFS, WINDOWS
                if (seed < save_max) {
                    gif.Close();
                    if (!headless) {
                        win.Close();
                    }
                }
            }
        }

        return;
    }


    /*

        DrawProjection()
            - two dimensional side view (a projection, not a slice) from z-camera POV

     */

    public static void DrawProjection(UIGrid vis, Parameters3D p, ExampleModel3d model, int time){
        // from z = 0 side
        System.out.println(p.sideLen);
        for (int xx = 0; xx < p.sideLen; xx++) {
            for (int yy = 0; yy < p.sideLen; yy++) {
                for (int zz = 0; zz < p.sideLen; zz++) {
                    Cell3D c=model.GetAgent(xx,yy,zz);
                    if(c!=null){
                        vis.SetPix(xx+2,yy+2,HSBColor(c.genome.h,c.genome.v,c.genome.s));
                        break;
                    } else {
                        vis.SetPix(xx+2, yy+2,Util.WHITE);
                    }
                }
            }
        }
//        vis.SetString(Integer.toString((int) time), true, BLACK, WHITE, 2);
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
