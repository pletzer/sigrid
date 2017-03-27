/**
 * Test command line argument parser
 */
#include <string>
#include <iostream>
#include <CmdLineArgParser.h>


int main(int argc, char** argv) {

    CmdLineArgParser prsr;
    prsr.set("--src_nj", 41, "Number of source grid latitudes");
    prsr.set("--delta_lat", 30.0, "Latitude pole displacement in deg.");
    prsr.set("--expr", "2*x + y", "Expression of x and y.");
    prsr.set("--name", std::string("foo"), "Name.");
    prsr.set("--vec", "0, 1, 2., -3.14, -4.e+06, -5.24e-45", "A vector.");
    prsr.parse(argc, argv);

    std::cout << "src_nj = " << prsr.get<int>("--src_nj") << '\n';
    std::cout << "delta_lat = " << prsr.get<double>("--delta_lat") << '\n';
    std::cout << "expr = " << prsr.get<const char*>("--expr") << '\n';
    std::cout << "name = " << prsr.get<std::string>("--name") << '\n';
    std::vector<double> vec = prsr.get<std::vector<double> >("--vec");
    std::cout << "vec = ";
    for (size_t i = 0; i < vec.size(); ++i) std::cout << vec[i] << ", ";
    std::cout << '\n';

    if (prsr.get<bool>("-h") || prsr.get<bool>("--help")) {
        prsr.help();
    }

    return 0;
}