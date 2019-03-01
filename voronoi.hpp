#pragma once

// point pattern analysis: compute Voronoi diagram for a set
// of seed points, and the IMTs of the Voronoi cells.
// optionally, with periodic boundary conditions
// 2019 Sebastian Kapfer <sebastian.kapfer@fau.de>
// 2019 Fabian Schaller <physik@fabian-schaller.de>
// FIXME this does not have tests yet

#include "papaya2.hpp"
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_policies_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Voronoi_diagram_2.h>

namespace papaya2 {
// put stuff in its own namespace as CGAL spams a lot
namespace voronoi_impl {

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using DT = CGAL::Delaunay_triangulation_2<K>;
using AT = CGAL::Delaunay_triangulation_adaptation_traits_2<DT>;
using AP = CGAL::Delaunay_triangulation_caching_degeneracy_removal_policy_2<DT>;
using VD = CGAL::Voronoi_diagram_2<DT, AT, AP>;

struct VoronoiDiagram : VD
{
    using point_t = K::Point_2;

    double boxLx = 1., boxLy = 1.;
    bool periodic = false;
    std::vector<point_t> seeds;

    void add_seed(vec_t seed) { seeds.push_back(point_t(seed[0], seed[1])); }

    // check if all the seeds are within the box.
    // this is required for the periodic version.
    // returns the first seed _not_ in the box in *offender if non-null.
    bool seeds_in_box(point_t *offender)
    {
        for (auto p : seeds) {
            if (!(0 <= p.x() && p.x() < boxLx && 0 <= p.y() && p.y() < boxLy)) {
                if (offender)
                    *offender = p;
                return false;
            }
        }
        return true;
    }

    void make_voronoi()
    {
        this->clear();

        if (periodic) {
            for (auto p : seeds)
                this->insert(p);
            add_padding_seeds(1.);
        } else {
            for (auto p : seeds)
                this->insert(p);
        }

        if (!this->is_valid())
            throw std::runtime_error(
                "CGAL thinks the Voronoi diagram is invalid.");
    }

    MinkowskiAccumulator minkval_for_cell(VD::Face_iterator fit) const
    {
        // circulate counterclockwise over the boundary
        auto heh = fit->ccb();
        auto const start = heh;

        std::vector<vec_t> vertices;
        vertices.reserve(10);
        do {
            auto p = heh->source()->point();
            vertices.push_back(point_to_vec(p));
        } while (++heh != start);

        return imt_polygon(vertices);
    }

  private:
    // FIXME make this adaptive
    void add_padding_seeds(double to)
    {
        double bbx_low = -boxLx * to;
        double bbx_high = (1 + to) * boxLx;
        double bby_low = -boxLy * to;
        double bby_high = (1 + to) * boxLy;
        for (int i = -1; i <= 1; ++i)
            for (int j = -1; j <= 1; ++j) {
                if (i == 0 && j == 0)
                    continue;

                for (auto p : seeds) {
                    point_t pp = {p.x() + i * boxLx, p.y() + j * boxLy};
                    if (pp.x() < bbx_low)
                        continue;
                    if (pp.x() >= bbx_high)
                        continue;
                    if (pp.y() < bby_low)
                        continue;
                    if (pp.y() >= bby_high)
                        continue;
                    this->insert(pp);
                }
            }
    }

    static vec_t point_to_vec(const point_t &p) { return {p.x(), p.y()}; }
};

} // namespace voronoi_impl

struct VoronoiDiagram : voronoi_impl::VoronoiDiagram
{};

} // namespace papaya2
