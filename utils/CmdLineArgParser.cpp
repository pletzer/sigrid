
#include <CmdLineArgParser.h>

// standard includes
#include <cstdlib>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <sstream>

CmdLineArgParser::CmdLineArgParser() {
    // always activate -h
    this->set("-h", false, "Print help.");
    this->set("--help", false, "Print help.");
    this->footer = "\nReport bugs to alexander@gokliya.net\n";
}
  
CmdLineArgParser::~CmdLineArgParser() {}

void
CmdLineArgParser::addFootnote(const std::string& note) {
    this->footer = note + "\n" + this->footer;
}

void 
CmdLineArgParser::set(const std::string& name, 
                      double defaultVal, 
                      const std::string& help) {
    this->doubleArg[name] = defaultVal;
    this->doubleArgHelp[name] = help;
}
  
void 
CmdLineArgParser::set(const std::string& name, 
                      int defaultVal, 
                      const std::string& help) {
    this->intArg[name] = defaultVal;
    this->intArgHelp[name] = help;
}
  
void 
CmdLineArgParser::set(const std::string& name, 
                      const std::string& defaultVal, 
                      const std::string& help) {
    this->stringArg[name] = defaultVal;
    this->stringArgHelp[name] = help;
}

void 
CmdLineArgParser::set(const std::string& name, 
                      const char* defaultVal, 
                      const std::string& help) {
    this->stringArg[name] = std::string(defaultVal);
    this->stringArgHelp[name] = help;
}

void 
CmdLineArgParser::set(const std::string& name, 
                      bool defaultVal, 
                      const std::string& help) {
    this->boolArg[name] = defaultVal;
    this->boolArgHelp[name] = help;
}

bool
CmdLineArgParser::parse(int argc, char *argv[]) {
    bool isArgVal;
    std::string nm;
    this->execName = argv[0];

    // check if all the options supplied are valid
    for (int a = 1; a < argc; ++a) {
        std::string arg(argv[a]);
        // is arg the name of an option (or a value)?
        bool isOptionName = (arg[0] == '-') 
            && (arg.size() >= 2) && !isdigit(arg[1]);
        if (isOptionName) {
            bool isValidOptName = 
               doubleArg.find(arg) != doubleArg.end()
            || intArg.find(arg) != intArg.end()
            || stringArg.find(arg) != stringArg.end()
            || boolArg.find(arg) != boolArg.end();
            if (!isValidOptName) {
                std::cout << arg << " is not a valid option.\n";
                return false;
            }
        } 
    }
  
    for (std::map<std::string, double>::iterator 
         i = this->doubleArg.begin(); 
         i != this->doubleArg.end(); ++i) {
        std::string optKey = i->first;
    
        isArgVal = false;
        nm = "";
        for (int a = 1; a < argc; ++a) {
            std::string arg(argv[a]);
            if (isArgVal) {
                i->second = atof(arg.c_str());
                isArgVal = false;
            }
            if (optKey == arg) {
                isArgVal = true;
                nm = optKey;
            }
        }
    }
    for (std::map<std::string, int>::iterator 
         i = this->intArg.begin(); 
         i != this->intArg.end(); ++i) {
       std::string optKey = i->first;
    
       isArgVal = false;
       nm = "";
       for (int a = 1; a < argc; ++a) {
           std::string arg(argv[a]);
           if (isArgVal) {
               i->second = atoi(arg.c_str());
               isArgVal = false;
           }
           if (optKey == arg) {
               isArgVal = true;
               nm = optKey;
           }
        }
    }
    for (std::map<std::string, std::string>::iterator 
         i = this->stringArg.begin(); 
         i != this->stringArg.end(); ++i) {
        std::string optKey = i->first;
    
        isArgVal = false;
        nm = "";
        for (int a = 1; a < argc; ++a) {
            std::string arg(argv[a]);
            if (isArgVal) {
                i->second = arg;
                isArgVal = false;
            }
            if (optKey == arg) {
                isArgVal = true;
                nm = optKey;
            }
        }
    }
    for (std::map<std::string, bool>::iterator 
         i = this->boolArg.begin(); 
         i != this->boolArg.end(); ++i) {
        std::string optKey = i->first;
        for (int a = 1; a < argc; ++a) {
            std::string arg(argv[a]);
            if (arg == i->first) {
                // flip the value
                i->second = !i->second;
            }
        }
    }

    return true;
}

void
CmdLineArgParser::help() const {
    std::cout << this->purpose << std::endl;
    std::cout << this->execName << " [options]\n";
    std::cout << "Usage:\n";
    for (std::map<std::string, double>::const_iterator 
         i = this->doubleArg.begin(); 
         i != this->doubleArg.end(); ++i) {
        std::string hlp = this->doubleArgHelp.find(i->first)->second;
        std::cout << "\t" << i->first << " <double#> " << hlp 
              << " (" << i->second << ")\n";
    }
    for (std::map<std::string, int>::const_iterator 
         i = this->intArg.begin(); 
         i != this->intArg.end(); ++i) {
        std::string hlp = this->intArgHelp.find(i->first)->second;
        std::cout << "\t" << i->first << " <int#> " << hlp 
              << " (" << i->second << ")\n";
    }
    for (std::map<std::string, std::string>::const_iterator 
         i = this->stringArg.begin(); 
         i != this->stringArg.end(); ++i) {
        std::string hlp = this->stringArgHelp.find(i->first)->second;
        std::cout << "\t" << i->first << " <string> " << hlp 
              << " (" << i->second << ")\n";
    }
    for (std::map<std::string, bool>::const_iterator 
         i = this->boolArg.begin(); 
         i != this->boolArg.end(); ++i) {
        std::string hlp = this->boolArgHelp.find(i->first)->second;
        std::cout << "\t" << i->first << ' ' << hlp 
              << " (" << i->second << ")\n";
    }
    std::cout << this->footer << std::endl;
}

template <class T>
T 
CmdLineArgParser::get(const std::string& name) const {
    std::ostringstream error;
    error << 
      "CmdLineArgParser::get: ERROR. Invalid type, use get<TYPE>(...)\n";
    throw error.str();
}

template <>
double 
CmdLineArgParser::get<double>(const std::string& name) const {
    double res = -std::numeric_limits<double>::max();
    std::map<std::string, double>::const_iterator i = this->doubleArg.find(name);
    if (i != this->doubleArg.end()) {
        res = i->second;
    }
    return res;
}

template <>
int 
CmdLineArgParser::get<int>(const std::string& name) const {
    int res = -std::numeric_limits<int>::max();
    std::map<std::string, int>::const_iterator i = this->intArg.find(name);
    if (i != this->intArg.end()) {
        res = i->second;
    }
    return res;
}

template <>
std::string
CmdLineArgParser::get<std::string>(const std::string& name) const {
    std::string res = "";
    std::map<std::string, std::string>::const_iterator 
      i = this->stringArg.find(name);
    if (i != this->stringArg.end()) {
        res = i->second;
    }
    return res;
}

template <>
const char*
CmdLineArgParser::get<const char*>(const std::string& name) const {
    std::map<std::string, std::string>::const_iterator 
      i = this->stringArg.find(name);
    if (i != this->stringArg.end()) {
        return i->second.c_str();
    }
    return "";
}

template <>
bool
CmdLineArgParser::get<bool>(const std::string& name) const {
    bool res = false;
    std::map<std::string, bool>::const_iterator 
    i = this->boolArg.find(name);
    if (i != this->boolArg.end()) {
        res = i->second;
    }
    return res;
}

template <>
std::vector<double>
CmdLineArgParser::get<std::vector<double> >(const std::string& name) const {
    // construct a double vector from a string such as "1.0, 2.e-2, 3.e33"
    std::vector<double> res;
    std::map<std::string, std::string>::const_iterator 
      i = this->stringArg.find(name);
    if (i != this->stringArg.end() && i->second.size() > 0) {
        std::string str = i->second;
        size_t pos0 = 0;
        size_t pos1 = str.find(',', pos0);
        while (pos1 < str.size()){
            res.push_back(std::atof(str.substr(pos0, pos1 - pos0).c_str()));
            pos0 = pos1 + 1;
            pos1 = str.find(',', pos0);
        }
        res.push_back(std::atof(str.substr(pos0, pos1 - pos0).c_str()));
    }
    return res;
}

template <>
std::vector<int>
CmdLineArgParser::get<std::vector<int> >(const std::string& name) const {
    // construct an int vector from a string such as "1, -2, 3"
    std::vector<int> res;
    std::map<std::string, std::string>::const_iterator 
      i = this->stringArg.find(name);
    if (i != this->stringArg.end() && i->second.size() > 0) {
        std::string str = i->second;
        size_t pos0 = 0;
        size_t pos1 = str.find(',', pos0);
        while (pos1 < str.size()){
            res.push_back(std::atoi(str.substr(pos0, pos1 - pos0).c_str()));
            pos0 = pos1 + 1;
            pos1 = str.find(',', pos0);
        }
        res.push_back(std::atof(str.substr(pos0, pos1 - pos0).c_str()));
    }
    return res;
}