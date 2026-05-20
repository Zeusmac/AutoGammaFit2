#ifndef OPTIONALDEPS_H
#define OPTIONALDEPS_H

#include <cstdlib>
#include <string>

// Runtime detection of optional dependencies.
// AutoGammaFit requires only ROOT. Everything here enables bonus features.
struct OptionalDeps {
    bool python3    = false;   // python3 in PATH
    bool plotly     = false;   // python3 -c 'import plotly'
    bool matplotlib = false;   // python3 -c 'import matplotlib'
    bool gnuplot    = false;   // gnuplot in PATH
    bool xdgOpen    = false;   // xdg-open (open files/URLs in browser)

    static OptionalDeps Probe() {
        OptionalDeps d;
        d.xdgOpen    = (system("which xdg-open  > /dev/null 2>&1") == 0);
        d.gnuplot    = (system("which gnuplot    > /dev/null 2>&1") == 0);
        d.python3    = (system("python3 --version > /dev/null 2>&1") == 0);
        if (d.python3) {
            d.plotly     = (system("python3 -c 'import plotly'     > /dev/null 2>&1") == 0);
            d.matplotlib = (system("python3 -c 'import matplotlib' > /dev/null 2>&1") == 0);
        }
        return d;
    }

    std::string Summary() const {
        std::string s;
        if (python3)    s += " python3";
        if (plotly)     s += " plotly";
        if (matplotlib) s += " matplotlib";
        if (gnuplot)    s += " gnuplot";
        if (xdgOpen)    s += " xdg-open";
        return s.empty() ? "none" : s.substr(1);
    }
};

#endif // OPTIONALDEPS_H
