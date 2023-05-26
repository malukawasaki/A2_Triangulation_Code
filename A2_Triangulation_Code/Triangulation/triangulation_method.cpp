/**
 * Copyright (C) 2015 by Liangliang Nan (liangliang.nan@gmail.com)
 * https://3d.bk.tudelft.nl/liangliang/
 *
 * This file is part of Easy3D. If it is useful in your research/work,
 * I would be grateful if you show your appreciation by citing it:
 * ------------------------------------------------------------------
 *      Liangliang Nan.
 *      Easy3D: a lightweight, easy-to-use, and efficient C++
 *      library for processing and rendering 3D data. 2018.
 * ------------------------------------------------------------------
 * Easy3D is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License Version 3
 * as published by the Free Software Foundation.
 *
 * Easy3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "triangulation.h"
#include "matrix_algo.h"
#include <easy3d/optimizer/optimizer_lm.h>
#include <numeric>


using namespace easy3d;


/**
 * TODO: Finish this function for reconstructing 3D geometry from corresponding image points.
 * @return True on success, otherwise false. On success, the reconstructed 3D points must be written to 'points_3d'
 *      and the recovered relative pose must be written to R and t.
 */
bool Triangulation::triangulation(
        double fx, double fy,     /// input: the focal lengths (same for both cameras)
        double cx, double cy,     /// input: the principal point (same for both cameras)
        const std::vector<Vector2D> &points_0,  /// input: 2D image points in the 1st image.
        const std::vector<Vector2D> &points_1,  /// input: 2D image points in the 2nd image.
        std::vector<Vector3D> &points_3d,       /// output: reconstructed 3D points
        Matrix33 &R,   /// output: 3 by 3 matrix, which is the recovered rotation of the 2nd camera
        Vector3D &t    /// output: 3D vector, which is the recovered translation of the 2nd camera
) const {
    /// NOTE: there might be multiple workflows for reconstructing 3D geometry from corresponding image points.
    ///       This assignment uses the commonly used one explained in our lecture.
    ///       It is advised to define a function for the sub-tasks. This way you have a clean and well-structured
    ///       implementation, which also makes testing and debugging easier. You can put your other functions above
    ///       triangulation(), or put them in one or multiple separate files.

    std::cout << "\nTODO: I am going to implement the triangulation() function in the following file:" << std::endl
              << "\t    - triangulation_method.cpp\n\n";

    std::cout << "[Liangliang]:\n"
                 "\tFeel free to use any provided data structures and functions. For your convenience, the\n"
                 "\tfollowing three files implement basic linear algebra data structures and operations:\n"
                 "\t    - Triangulation/matrix.h  Matrices of arbitrary dimensions and related functions.\n"
                 "\t    - Triangulation/vector.h  Vectors of arbitrary dimensions and related functions.\n"
                 "\t    - Triangulation/matrix_algo.h  Determinant, inverse, SVD, linear least-squares...\n"
                 "\tPlease refer to the above files for a complete list of useful functions and their usage.\n\n"
                 "\tIf you choose to implement the non-linear method for triangulation (optional task). Please\n"
                 "\trefer to 'Tutorial_NonlinearLeastSquares/main.cpp' for an example and some explanations.\n\n"
                 "\tIn your final submission, please\n"
                 "\t    - delete ALL unrelated test or debug code and avoid unnecessary output.\n"
                 "\t    - include all the source code (and please do NOT modify the structure of the directories).\n"
                 "\t    - do NOT include the 'build' directory (which contains the intermediate files in a build step).\n"
                 "\t    - make sure your code compiles and can reproduce your results without ANY modification.\n\n"
              << std::flush;

    /// Below are a few examples showing some useful data structures and APIs.

    /// define a 2D vector/point
    Vector2D b(1.1, 2.2);

    /// define a 3D vector/point
    Vector3D a(1.1, 2.2, 3.3);

    /// get the Cartesian coordinates of a (a is treated as Homogeneous coordinates)
    Vector2D p = a.cartesian();

    /// get the Homogeneous coordinates of p
    Vector3D q = p.homogeneous();

    /// define a 3 by 3 matrix (and all elements initialized to 0.0)
    Matrix33 A;

    /// define and initialize a 3 by 4 matrix
    Matrix34 M(1.1, 2.2, 3.3, 0,
               0, 2.2, 3.3, 1,
               0, 0, 1, 1);

    /// set first row by a vector
    M.set_row(0, Vector4D(1.1, 2.2, 3.3, 4.4));

    /// set second column by a vector
    M.set_column(1, Vector3D(5.5, 5.5, 5.5));

    /// define a 15 by 9 matrix (and all elements initialized to 0.0)
    Matrix W(15, 9, 0.0);
    /// set the first row by a 9-dimensional vector
    W.set_row(0, {0, 1, 2, 3, 4, 5, 6, 7, 8}); // {....} is equivalent to a std::vector<double>

    /// get the number of rows.
    int num_rows = W.rows();

    /// get the number of columns.
    int num_cols = W.cols();

    /// get the the element at row 1 and column 2
    double value = W(1, 2);

    /// get the last column of a matrix
    Vector last_column = W.get_column(W.cols() - 1);

    /// define a 3 by 3 identity matrix
    Matrix33 I = Matrix::identity(3, 3, 1.0);

    /// matrix-vector product
    Vector3D v = M * Vector4D(1, 2, 3, 4); // M is 3 by 4

    ///For more functions of Matrix and Vector, please refer to 'matrix.h' and 'vector.h'

    // TODO: delete all above example code in your final submission

    //--------------------------------------------------------------------------------------------------------------
    // implementation starts ...

    /// Check if the input is valid (i.e. sizes of points must match and be greater or equal to 8)
    if (points_0.size() != points_1.size() || points_0.size() < 8) {
        std::cout << "Sizes of points do not match or are smaller than 8. This operation is not possible" << std::endl;
        return false;
    } else {
        std::cout << "Sizes of points match and are greater or equal to 8. This operation is possible"
                  << std::endl;

        // TODO: Estimate relative pose of two views. This can be done by solving the following steps:
        //      - estimate the fundamental matrix F;
        //      - compute the essential matrix E;
        //      - recover rotation R and t.

        // TODO: Estimate the fundamental matrix F;
        // Normalization
        double sumx0 = 0;
        double sumy0 = 0;
        double sumx1 = 0;
        double sumy1 = 0;
        for (int i = 0; i < points_0.size(); i++) {
            sumx0 = sumx0 + points_0[i].x();
            sumy0 = sumy0 + points_0[i].y();
            sumx1 = sumx1 + points_1[i].x();
            sumy1 = sumy1 + points_1[i].y();}

        double tx0 = sumx0 / points_0.size();
        double ty0 = sumy0 / points_0.size();
        double tx1 = sumx1 / points_1.size();
        double ty1 = sumy1 / points_1.size();
        Vector2D s0;
        Vector2D s1;

        for (int i = 0; i < points_0.size(); i++) {
            Vector2D mc0[i];
            Vector2D mc1[i];
            Vector2D dc0[i];
            Vector2D dc1[i];
            mc0[i].x() = points_0[i].x() - tx0;
            mc0[i].y() = points_0[i].y() - ty0;
            mc1[i].x() = points_1[i].x() - tx1;
            mc1[i].y() = points_1[i].y() - ty1;
            dc0[i].x() += sqrt(pow(mc0[i].x(),2));
            dc0[i].y() += sqrt(pow(mc0[i].y(),2));
            dc1[i].x() += sqrt(pow(mc1[i].x(),2));
            dc1[i].y() += sqrt(pow(mc1[i].y(),2));
            s0[i] = sqrt(2)/dc0[i];
            //TODO: Check if the one bellow is a way calculate s
            s1[i] = dot(1.0/sqrt(2),dc1[i]);
        }

        // TODO: Fix s
        // TODO: Do we need two matrices? (i.e. T0 and T1 to use tx0 and tx1?)
        Matrix33 T0(s0, 0, -s0*tx0,
                    0, s0, -s0*ty0,
                    0, 0, 1);

        Matrix33 T1(s1, 0, -s1*tx1,
                    0, s1, -s1*ty1,
                    0, 0, 1);

        //F = F'*T0*T1?

        // TODO: Reconstruct 3D points. The main task is
        //      - triangulate a pair of image points (i.e., compute the 3D coordinates for each corresponding point pair)

        // TODO: Don't forget to
        //          - write your recovered 3D points into 'points_3d' (so the viewer can visualize the 3D points for you);
        //          - write the recovered relative pose into R and t (the view will be updated as seen from the 2nd camera,
        //            which can help you check if R and t are correct).
        //       You must return either 'true' or 'false' to indicate whether the triangulation was successful (so the
        //       viewer will be notified to visualize the 3D points and update the view).
        //       There are a few cases you should return 'false' instead, for example:
        //          - function not implemented yet;
        //          - input not valid (e.g., not enough points, point numbers don't match);
        //          - encountered failure in any step.

        return points_3d.size() > 0;
    }
}