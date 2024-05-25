﻿#include <iostream>
#include <fstream>
#include <iomanip>
#include "saspoint5_distribution.hpp"

void plot_pdf(std::string filepath) {
    ofstream ofs(filepath);

    ofs << "x,pdf" << std::endl;
    ofs << std::scientific << std::setprecision(16);

    for (double x = -6; x <= 64; x += 1. / 1024) {
        ofs << x << ',' << saspoint5_pdf(x) << std::endl;
    }

    ofs.close();
}

void plot_pdf_limit(std::string filepath) {
    ofstream ofs(filepath);

    ofs << "x,pdf" << std::endl;
    ofs << std::scientific << std::setprecision(16);

    for (double x0 = 64; x0 <= ldexp(1, 64); x0 *= 2) {
        for (double x = x0; x < x0 * 2; x += x0 / 256) {
            ofs << x << ',' << saspoint5_pdf(x) << std::endl;
        }
    }

    ofs.close();
}

void plot_cdf(std::string filepath) {
    ofstream ofs(filepath);

    ofs << "x,cdf,ccdf" << std::endl;
    ofs << std::scientific << std::setprecision(16);

    for (double x = -6; x <= 64; x += 1. / 1024) {
        ofs << x << ',' << saspoint5_cdf(x) << ',' << saspoint5_cdf(x, true) << std::endl;
    }

    ofs.close();
}

void plot_cdf_limit(std::string filepath) {
    ofstream ofs(filepath);

    ofs << "x,ccdf" << std::endl;
    ofs << std::scientific << std::setprecision(16);

    for (double x0 = 64; x0 <= ldexp(1, 64); x0 *= 2) {
        for (double x = x0; x < x0 * 2; x += x0 / 256) {
            ofs << x << ',' << saspoint5_cdf(x, true) << std::endl;
        }
    }

    ofs.close();
}

void plot_quantile(std::string filepath) {
    ofstream ofs(filepath);

    ofs << "x,quantile" << std::endl;
    ofs << std::scientific << std::setprecision(16);

    for (double x = 1. / 8192; x < 1; x += 1. / 8192) {
        ofs << x << ',' << saspoint5_quantile(x) << std::endl;
    }

    ofs.close();
}

void plot_quantilelower_limit(std::string filepath) {
    ofstream ofs(filepath);

    ofs << "x,quantile" << std::endl;
    ofs << std::scientific << std::setprecision(16);

    for (double x = 1. / 8192; x > ldexp(1, -1000); x /= 2) {
        ofs << x << ',' << saspoint5_quantile(x) << std::endl;
    }

    ofs.close();
}

void plot_quantileupper_limit(std::string filepath) {
    ofstream ofs(filepath);

    ofs << "x,cquantile" << std::endl;
    ofs << std::scientific << std::setprecision(16);

    for (double x0 = 1. / 8192; x0 > ldexp(1, -128); x0 /= 2) {
        for (double x = x0; x > x0 / 2; x -= x0 / 256) {
            ofs << x << ',' << saspoint5_quantile(x, true) << std::endl;
        }
    }

    ofs.close();
}

int main() {
    plot_pdf("../results/saspoint5_pdf_cpp.csv");
    plot_pdf_limit("../results/saspoint5_pdf_limit_cpp.csv");
    plot_cdf("../results/saspoint5_cdf_cpp.csv");
    plot_cdf_limit("../results/saspoint5_cdf_limit_cpp.csv");
    plot_quantile("../results/saspoint5_quantile_cpp.csv");
    plot_quantilelower_limit("../results/saspoint5_quantilelower_limit_cpp.csv");
    plot_quantileupper_limit("../results/saspoint5_quantileupper_limit_cpp.csv");

    std::cout << "END" << std::endl;
}