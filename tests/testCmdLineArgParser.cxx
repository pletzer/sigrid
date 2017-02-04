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
    prsr.parse(argc, argv);

    std::cout << "src_nj = " << prsr.get<int>("--src_nj") << '\n';
    std::cout << "delta_lat = " << prsr.get<double>("--delta_lat") << '\n';
    std::cout << "expr = " << prsr.get<const char*>("--expr") << '\n';
    std::cout << "name = " << prsr.get<std::string>("--name") << '\n';

    if (prsr.get<bool>("-h") || prsr.get<bool>("--help")) {
        prsr.help();
    }

    return 0;
}