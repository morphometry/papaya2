#include "catch.hpp"
#include "tools.hpp"

using namespace papaya2;

TEST_CASE("fsq", "[tools]")
{
    CHECK(fsq(-2) == 4.);
    CHECK(std::isnan(fsq(NAN)));
}

TEST_CASE("space_separated and comma_separated works")
{
    std::vector<std::string> dummyvec = {"a", "b", "c"};
    SECTION("space_separated")
    {
        std::ostringstream os;
        os << space_separated(dummyvec);
        CHECK(os.str() == "a b c");
    }
    SECTION("comma_separated")
    {
        std::ostringstream os;
        os << comma_separated(dummyvec);
        CHECK(os.str() == "a, b, c");
    }
}

TEST_CASE("vector operations", "[vector]")
{
    SECTION("vector type promotion")
    {
        std::array<int, 2> int_array = {1, 3};
        std::array<double, 2> flt_array = {.5, 1.5};
        CHECK(.5 * int_array == flt_array);
        CHECK(int_array * .5 == flt_array);
        CHECK(int_array / 2. == flt_array);
    }

    SECTION("vector +=")
    {
        std::array<double, 2> flt_array = {1.5, 1.5};
        int const x[2] = {0, 2};
        flt_array += x;
        CHECK(flt_array[0] == 1.5);
        CHECK(flt_array[1] == 3.5);
    }
}

TEST_CASE("integer range container", "[range]")
{
    using intvec_t = std::vector<int>;

    SECTION("range-based for loop")
    {
        intvec_t v;
        for (int x : range(4, 6))
            v.push_back(x);
        CHECK(v == intvec_t({4, 5}));
    }

    SECTION("iterator interface")
    {
        intvec_t v;
        auto r = range(4);
        v.assign(r.begin(), r.end());
        CHECK(v == intvec_t({0, 1, 2, 3}));
    }
}
