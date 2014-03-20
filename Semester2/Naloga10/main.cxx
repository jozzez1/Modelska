#include <iostream>
#include <fstream>
#include <string>
#include <mgl2/data.h>
#include <mgl2/mgl.h>

int main (int argc, char ** argv)
{
    if (argc != 2)
    {
        std::cerr << "Wrong program usage! Try:" 	<< std::endl;
        std::cerr << argv[0] << " <filename>" 	    << std::endl;
        exit (EXIT_FAILURE);
    }

    // we get the filename
    std::string filename (argv[1]);

    // from the format we read N
    // fortmat is -N-md.txt ...
    size_t start    = filename.find_first_of ("-"),
           stop     = filename.find_last_of  ("-"),
           dot      = filename.find_last_of (".");

    // sp this the int
    size_t N = atoi ((filename.substr(start+1, stop-1)).c_str()) + 2;
    std::string title = "Poves, \\Delta m = " + filename.substr(stop+1,dot);
    title.replace (title.end()-4, title.end(), "");

    mglData x (N, N),
            y (N, N),
            z (N, N),
            X (13),
            Y (13);

    // let's also plot the outline of the weighted object
    X.SetVal (0.25, 0); Y.SetVal (0, 0);
    X.SetVal (0.25, 1); Y.SetVal (0.5, 1);
    X.SetVal (0, 2);    Y.SetVal (0.5, 2);
    X.SetVal (0, 3);    Y.SetVal (0.75, 3);
    X.SetVal (0.25, 4); Y.SetVal (0.75, 4);
    X.SetVal (0.25, 5); Y.SetVal (1, 5);
    X.SetVal (0.5, 6);  Y.SetVal (1, 6);
    X.SetVal (0.5, 7);  Y.SetVal (0.75, 7);
    X.SetVal (1, 8);    Y.SetVal (0.75, 8);
    X.SetVal (1, 9);    Y.SetVal (0.5, 9);
    X.SetVal (0.5, 10); Y.SetVal (0.5, 10);
    X.SetVal (0.5, 11); Y.SetVal (0, 11);
    X.SetVal (0.25, 12);Y.SetVal (0, 12);


    std::ifstream fin (filename.c_str(), std::ifstream::in);
    if (!fin)
    {
        std::cerr << "Error. File not found." << std::endl;
        exit (EXIT_FAILURE);
    }

    for (long i = 0; !fin.eof(); i++)
        fin >> x.a[i] >> y.a[i] >> z.a[i];
    fin.close ();

    // now we change the ending for output
    filename.replace (filename.end()-3, filename.end(), "png");

    mglGraph gr (0, 800, 800);
    gr.SetRanges (x,y);
    gr.SetRange ('c', z);
    gr.Colorbar ("wyqrRk_");
    gr.Axis("q");
    gr.Box ("q");

    gr.Cont (x, y, z, "t");
    gr.Dens (x, y, z, "wyqrRk");
    gr.Plot (X, Y, "W2i");
    gr.Title (title.c_str());

    gr.WritePNG (filename.c_str(), "Poves");
    gr.ShowImage ("feh");
    exit (EXIT_SUCCESS);
}

