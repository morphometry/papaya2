// point pattern analysis: compute Voronoi diagram for a set
// of seed points, and the IMTs of the Voronoi cells.
// optionally, with periodic boundary conditions
// 2019 Sebastian Kapfer <sebastian.kapfer@fau.de>
// 2019 Fabian Schaller <physik@fabian-schaller.de>
// FIXME this does not have tests yet
#include "common.hpp"
#ifndef HAVE_CGAL
#error ppanalysis requiresm the CGAL library to work.  (Set CGAL_SUPPORT = 1.)
#endif
#include <papaya2/voronoi.hpp>
using namespace papaya2;
using papaya2::VoronoiDiagram;
using point_t = VoronoiDiagram::point_t;

static std::vector<point_t> parse_plain_ascii_file(string infilename)
{
    std::vector<point_t> ret;
    std::ifstream ifs(infilename);
    double x, y;
    while (ifs >> x >> y)
        ret.push_back(point_t(x, y));
    return ret;
}

int main(int argc, const char **argv)
{
    string infilename, outfilename = "/dev/stdout";
    VoronoiDiagram vd;

    // process command-line arguments
    for (++argv; *argv; ++argv) {
        if (string(*argv) == "in") {
            infilename = read_arg<string>(argv++);
        } else if (string(*argv) == "out") {
            outfilename = read_arg<string>(argv++);
        } else if (string(*argv) == "boxL") {
            vd.boxLx = vd.boxLy = read_arg<double>(argv++);
            vd.periodic = true;
        } else if (string(*argv) == "-h") {
            std::cerr << "usage:\n (nonperiodic mode)\n\tppanalysis in <input "
                         "filename> out <output filename>\n";
            std::cerr << " (periodic mode)\n";
            std::cerr << "\tppanalysis in <input filename> out <output "
                         "filename> boxL <box side length>\n";
            std::cerr << "for help, please see "
                         "https://morphometry.org/software/papaya2/\n";
            return 0;
        } else {
            std::cerr << "illegal argument: " << *argv << "\n";
            return 1;
        }
    }

    if (infilename == "") {
        std::cerr << "input filename not set (use the 'in' argument)\n";
        return 1;
    }

    if (outfilename == "") {
        std::cerr << "output filename not set (use the 'out' argument)\n";
        return 1;
    }

    std::ofstream outfile(outfilename);
    outfile << "# IMTs of Voronoi cells of a point pattern ";
    if (vd.periodic)
        outfile << "with periodic boundary conditions\n# box size: " << vd.boxLx
                << " " << vd.boxLy;
    outfile << "\n# seedx seedy area perimeter q2 q3 q4 q5 q6\n";
    outfile.precision(8);
    if (!outfile) {
        std::cerr << "unable to open output file: " << outfilename << "\n";
        return 1;
    }

    vd.seeds = parse_plain_ascii_file(infilename);

    if (vd.periodic) {
        point_t offender;
        if (!vd.seeds_in_box(&offender)) {
            std::cerr << "There are points outside of the box: " << offender.x()
                      << ", " << offender.y() << "\n";
            return 1;
        }
    }

    vd.make_voronoi();

    int label = 0;
    for (auto fit = vd.faces_begin(); fit != vd.faces_end(); ++fit) {
        const point_t seed = fit->dual()->point();

        if (seed != vd.seeds.at(label)) {
            std::cerr << "CGAL decided to reorder our seed points.\n";
            return 1;
        }

        if (!fit->is_unbounded()) {
            auto minkval = vd.minkval_for_cell(fit);

            outfile << seed.x() << " " << seed.y() << " " << minkval.area()
                    << " " << minkval.perimeter();
            for (auto s : {2, 3, 4, 5, 6}) {
                outfile << " " << minkval.msm(s);
            }
            outfile << "\n";
        } else {
            outfile << seed.x() << " " << seed.y() << " unbounded\n";
        }

        ++label;
        if (label == (int)vd.seeds.size())
            break;
    }

    return 0;
}
