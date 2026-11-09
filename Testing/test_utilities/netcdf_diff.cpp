//------------------------------------------------------------------------------
//  The @header2, @begin_table, @item3 and @end_table commands are custom
//  defined commands in Doxygen.in. They are defined under ALIASES. For the page
//  created here, the 80 column limit is exceeded. Arguments of aliases are
//  separated by ','. If you intended ',' to be a string you must use an escaped
//  comma '\,'.
//
///  @page netcdf_diff_cl_parsing_sec Command Line Arguments
///
///  @tableofcontents
///
///  @section netcdf_diff_cl_parsing_intro Introduction
///  This contains a description of the command line arguments. All arguments
///  take the form of
///
///  @fixed_width{-arg=value}
///
///  @section netcdf_diff_cl_parsing_arg_sec Command Line Arguments
///  @header2{Argument, Takes Value, Discription}
///  @begin_table
///     @item3{@fixed_width{-h},            N, Displays the help text and exits the program.}
///     @item3{@fixed_width{-netcdf_file1}, Y, First netcdf file name.}
///     @item3{@fixed_width{-netcdf_file2}, Y, Second netcdf file name.}
///     @item3{@fixed_width{-quantity},     Y, Netcdf quantity to check.}
///     @item3{@fixed_width{-tol},          Y, Tolarance value.}
///  @end_table
///
///  @section netcdf_diff_cl_pasring_prog_ref_sec Programmers Reference
///  Reference material for the coding to implement command line parsing is found
///  in the @ref netcdf_test::commandline_parser class.
//------------------------------------------------------------------------------
//******************************************************************************
///  @file netcdf_diff.cpp
///  @brief Utility to check the difference between two netcdf files.
//******************************************************************************

#include <cmath>
#include <iostream>
#include <vector>
#include <netcdf.h>
#include <algorithm>

#include "netcdf_test_commandline_parser.hpp"

//------------------------------------------------------------------------------
///  @brief Load a quantity from a netcdf file.
///
///  Quantites are loaded into a flat vector for rapid comparison.
///
///  @param[in] netcdf_file Netcdf file name.
///  @param[in] name        Name of the netcdf file quantity.
//------------------------------------------------------------------------------
std::vector<double> netcdf_quantity(const std::string netcdf_file,
                                    const std::string name) {
    int ncid;
    nc_open(netcdf_file.c_str(), NC_NOWRITE, &ncid);

    int varid;
    nc_inq_varid(ncid, name.c_str(), &varid);

    int ndims;
    nc_inq_varndims(ncid, varid, &ndims);

    std::vector<int> dimids(ndims);
    nc_inq_vardimid(ncid, varid, dimids.data());

    size_t total_length = 1;
    for (int dimid: dimids) {
        size_t dim_length;
        nc_inq_dimlen(ncid, dimid, &dim_length);

        total_length *= dim_length;
    }

    std::vector<double> buffer(total_length);
    nc_get_var(ncid, varid, buffer.data());

    nc_close(ncid);

    return buffer;
}

//------------------------------------------------------------------------------
///  @brief Main test program.
///
///  @param[in] argc Number of commandline arguments.
///  @param[in] argv Array of arguments strings.
//------------------------------------------------------------------------------
int main(int argc, const char * argv[]) {
    const netcdf_test::commandline_parser args(argc, argv, []() -> void {
//                   "                                        ''                                      "
        std::cout << "                                                                                " << std::endl;
        std::cout << "                                    NETCDF DIFF                                 " << std::endl;
        std::cout << "                                                                                " << std::endl;
        std::cout << "Usage: netcdf_diff [-arg][=option] ...                                          " << std::endl;
        std::cout << "                                                                                " << std::endl;
        std::cout << "Options:                                                                        " << std::endl;
        std::cout << "All options are displayes as [arg][takesoption][Discription]                    " << std::endl;
        std::cout << "  -h            N Display this information.                                     " << std::endl;
        std::cout << "                                                                                " << std::endl;
        std::cout << "  -netcdf_file1 Y First netcdf file name.                                       " << std::endl;
        std::cout << "                                                                                " << std::endl;
        std::cout << "  -netcdf_file2 Y Second netcdf file name.                                      " << std::endl;
        std::cout << "                                                                                " << std::endl;
        std::cout << "  -quantity     Y Netcdf quantity to check.                                     " << std::endl;
        std::cout << "                                                                                " << std::endl;
        std::cout << "  -tol          Y Tolarance value.                                              " << std::endl;
        std::cout << "                                                                                " << std::endl;
        std::cout << std::endl;

        exit(1);
    });

    std::string quantity = args.get<std::string> ("-quantity");

    std::vector<double> q1 =
        netcdf_quantity(args.get<std::string> ("-netcdf_file1"), quantity);
    std::vector<double> q2 =
        netcdf_quantity(args.get<std::string> ("-netcdf_file2"), quantity);

    const double tolarance = args.get<double> ("-tol");

    bool pass = q1.size() == q2.size();
    if (!pass) {
        std::cout << "Quantity " << quantity
                  << " has unequal lengths." << std::endl;
        std::cout << q1.size() << " " << q2.size() << std::endl;
        exit(1);
    }

    for (size_t i = 0, e = q1.size(); i < e; i++) {
        pass = pass && std::abs(q1[i] - q2[i]) < tolarance;
        if (!pass) {
            std::cout << std::abs(q1[i] - q2[i]) << std::endl;
        }
    }

    if (!pass) {
        std::cout << "Quantity " << quantity
                  << " has unequal values." << std::endl;
        exit(1);
    }

    exit(0);
}
