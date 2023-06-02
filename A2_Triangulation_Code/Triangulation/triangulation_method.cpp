/**
 * GEO10016 Assignment 1: Camera Calibration
 * Group 06: Maria Luisa Tarozzo Kawasaki (5620341), Simay Batum (5715598), Rianne Aalders (4593987)
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
#include <list>

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

    /*  /// define a 2D vector/point
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
      // Matrix W(15, 9, 0.0);
      /// set the first row by a 9-dimensional vector
      // W.set_row(0, {0, 1, 2, 3, 4, 5, 6, 7, 8}); // {....} is equivalent to a std::vector<double>

      /// get the number of rows.
     // int num_rows = W.rows();

      /// get the number of columns.
      // int num_cols = W.cols();

      /// get the the element at row 1 and column 2
     // double value = W(1, 2);

      /// get the last column of a matrix
      // Vector last_column = W.get_column(W.cols() - 1);

      /// define a 3 by 3 identity matrix
      Matrix33 I = Matrix::identity(3, 3, 1.0);

      /// matrix-vector product
      Vector3D v = M * Vector4D(1, 2, 3, 4); // M is 3 by 4

      ///For more functions of Matrix and Vector, please refer to 'matrix.h' and 'vector.h'*/

    // TODO: delete all above example code in your final submission

    //--------------------------------------------------------------------------------------------------------------

    /// Check if the input is valid (i.e. sizes of points must match and be greater or equal to 8)
    if (points_0.size() != points_1.size() || points_0.size() < 8) {
        std::cout << "Sizes of points do not match or are smaller than 8. This operation is not possible" << std::endl;
        return false;
    } else {
        std::cout << "Sizes of points match and are greater or equal to 8. This operation is possible"
                  << std::endl;

        // STEP 1. Estimate the fundamental matrix F.
        //Calculate the centers (corresponding pixel centers in each image separately).
        Vector2D center0(0,0);
        Vector2D center1(0,0);

        for (int i = 0; i < points_0.size(); i++){
            center0+= points_0[i];
            center1 += points_1[i];
        }
        center0 /= points_0.size();
        center1 /= points_1.size();

        //Compute mean distance to center (in each image separately).
        double avg_distance0 = 0;
        double avg_distance1 = 0;
        for (int i = 0; i < points_0.size(); i++){
            //Adds the value of points1[i] to the "avg_distance1" variable.
            avg_distance0 += (points_0[i] - center0).norm();
            avg_distance1 += (points_1[i] - center1).norm();
        }

        avg_distance0 /= points_0.size();
        avg_distance1 /= points_1.size();

        /*std::cout << "ADis0 = " << avg_distance0 << std::endl;
        std::cout << "Adis1 = " << avg_distance1 << std::endl;*/

        //Compute scaling factor using average distance.
        double scaling_factor0 = sqrt(2)/avg_distance0;
        double scaling_factor1 = sqrt(2)/avg_distance1;

        //Compute initial Fundamental matrix with norm points.
        std::vector<Vector2D> points_0_normalized(points_0.size());
        std::vector<Vector2D> points_1_normalized(points_1.size());
        for (int i = 0; i < points_0.size(); i++){
            points_0_normalized[i].x() = (points_0[i].x() - center0.x()) * scaling_factor0;
            points_0_normalized[i].y() = (points_0[i].y() - center0.y()) * scaling_factor0; // I changed the x's to y's in this line
            points_1_normalized[i].x() = (points_1[i].x() - center1.x()) * scaling_factor1;
            points_1_normalized[i].y() = (points_1[i].y() - center1.y()) * scaling_factor1; // I changed the x's to y's in this line
        }

        int size = int(points_0.size());
        Matrix W (size,9);
        //W = [xi'xi, xi'yi, xi', yi'xi, yi'yi, yi', xi, yi, 1]
        for (int i = 0; i < points_0.size(); i++){
            W.set(i,0, points_0_normalized[i].x() * points_1_normalized[i].x());
            W.set(i,1, points_0_normalized[i].y() * points_1_normalized[i].x());
            W.set(i,2, points_1_normalized[i].x());
            W.set(i,3, points_0_normalized[i].x() * points_1_normalized[i].y());
            W.set(i,4, points_0_normalized[i].y() * points_1_normalized[i].y());
            W.set(i,5, points_1_normalized[i].y());
            W.set(i,6, points_0_normalized[i].x());
            W.set(i,7, points_0_normalized[i].y());
            W.set(i,8, 1);
        }

        //Extract the right singular vector corresponding to the smallest singular value. Reshape it into a 3x3 matrix F_hat.
        int n_cols = W.cols();
        int n_rows = W.rows();
        Matrix U(n_rows, n_rows);
        Matrix V(n_cols,n_cols);
        Matrix D(n_rows,n_cols);
        Matrix33 F_hat;

        svd_decompose(W, U, D, V);

        //Populating F_hat with the last column of V.
        int col_index = 8;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                F_hat(i, j) = V(i * 3 + j, col_index);
            }
        }
        //std::cout << "F_hat = " << F_hat << std::endl;

        Matrix33 A, B, C;
        svd_decompose(F_hat, A, B, C);

        // Rank 2 enforcement:
        // We can enforce the rank-2 constraint by setting the smallest singular value to 0.
        B.set(2,2,0);
        //std::cout << "B = " << B << std::endl;

        Matrix33 Fq = A * B * C.transpose();
        //std::cout << "Fq = " << Fq << std::endl;

        // Calculate T0, and T1.
        Matrix33 T0(sqrt(2) / avg_distance0, 0, -sqrt(2) / avg_distance0 * center0.x(),
                    0, sqrt(2) / avg_distance0, -sqrt(2) / avg_distance0 * center0.y(),
                    0, 0, 1);
        Matrix33 T1(sqrt(2) / avg_distance1, 0, -sqrt(2) / avg_distance1 * center1.x(),
                    0, sqrt(2) / avg_distance1, -sqrt(2) / avg_distance1 * center1.y(),
                    0, 0, 1);

        // Denormalize Fq to F.
        Matrix33 F = T1.transpose() * Fq * T0;

        // STEP 2. Recover relative pose (i.e., R and t) from the fundamental matrix.
        // Calculate E.
        Matrix33 K (fx, 0, cx, // skew = 0 because the cross product is skew symmetric
                    0, fy, cy,
                    0,0,1);
        Matrix33 E = K.transpose() * F * K;
        //std::cout << "E = " << E << std::endl;
        Matrix33 X,Z;
        // Enforce that the diagonal of Ydiag = (1,1,0)
        Matrix33 Y(1,0,0,
                   0,1,0,
                   0,0,0);
        //std::cout << "Y = " << Y << std::endl;
        svd_decompose(E, X, Y, Z); // X=U, Y=D, Z=V
        //std::cout << "X = " << X << std::endl;
        //std::cout << "Z = " << Z << std::endl;

        // Calculate 4 Rt settings from E.
        Matrix33 E_W(0, -1, 0,
                     1, 0, 0,
                     0, 0, 1);
        Matrix33 E_Z(0, 1, 0,
                     -1, 0, 0,
                     0, 0, 0);

        Matrix33 tx = X * E_Z * X.transpose();
        std::cout << "tx = " << tx << std::endl;

        Matrix33 R1 = X * E_W *Z.transpose();
        Matrix33 R2 = X * E_W.transpose() *Z.transpose();
        std::cout << "R1 = " << R1 << std::endl;
        std::cout << "R2 = " << R2 << std::endl;

        R1 *= determinant(R1);
        R2 *= determinant(R2);
        std::cout << "determinant of R1 = " << determinant(R1) << std::endl;
        std::cout << "determinant of R2 = " << determinant(R2) << std::endl;

        Vector3D t1 = X.get_column(X.cols() - 1);
        Vector3D t2 = -X.get_column(X.cols() - 1);
        std::cout << "t1 = " << t1 << std::endl;
        std::cout << "t2 = " << t2 << std::endl;

        // Create matrices Rt with the 4 different combinations.
        Matrix34 Rt1(R1(0,0), R1(0,1),R1(0,2), t1.x(),
                     R1(1,0), R1(1,1),R1(1,2), t1.y(),
                     R1(2,0), R1(2,1),R1(2,2), t1.z());

        Matrix34 Rt2(R1(0,0), R1(0,1),R1(0,2), t2.x(),
                     R1(1,0), R1(1,1),R1(1,2), t2.y(),
                     R1(2,0), R1(2,1),R1(2,2), t2.z());

        Matrix34 Rt3(R2(0,0), R2(0,1),R2(0,2), t1.x(),
                     R2(1,0), R2(1,1),R2(1,2), t1.y(),
                     R2(2,0), R2(2,1),R2(2,2), t1.z());

        Matrix34 Rt4(R2(0,0), R2(0,1),R2(0,2), t2.x(),
                     R2(1,0), R2(1,1),R2(1,2), t2.y(),
                     R2(2,0), R2(2,1),R2(2,2), t2.z());


        // M0 refers to M = K[I0].
        Matrix34 M0_dummy(1, 0, 0, 0,
                          0, 1, 0, 0,
                          0, 0, 1, 0);

        Matrix34 M0;
        M0 = K * M0_dummy;

        // M1, M2, M3 and M4 refer to M' = K * Rt.
        Matrix34 M1 = K * Rt1;
        Matrix34 M2 = K * Rt2;
        Matrix34 M3 = K * Rt3;
        Matrix34 M4 = K * Rt4;

        // STEP 3. Create matrices A = [xm3^T - m1^T, ym3^T - m2^T, x'm'3^T - m'1^T, y'm'3^T - m'2^T].
        // Create 3D points using the different As and Ms.
        std::vector<Vector3D> P1_array, P2_array, P3_array, P4_array;
        Vector3D P1_3D, P2_3D, P3_3D, P4_3D;

        for (int i = 0; i < points_0.size(); i++) {
            Matrix44 A1, A2, A3, A4;
            A1.set_row(0, points_0[i].x() * M0.get_row(2) - M0.get_row(0));
            A1.set_row(1, points_0[i].y() * M0.get_row(2) - M0.get_row(1));
            A1.set_row(2, points_1[i].x() * M1.get_row(2) - M1.get_row(0));
            A1.set_row(3, points_1[i].y() * M1.get_row(2) - M1.get_row(1));

            A2.set_row(0, points_0[i].x() * M0.get_row(2) - M0.get_row(0));
            A2.set_row(1, points_0[i].y() * M0.get_row(2) - M0.get_row(1));
            A2.set_row(2, points_1[i].x() * M2.get_row(2) - M2.get_row(0));
            A2.set_row(3, points_1[i].y() * M2.get_row(2) - M2.get_row(1));

            A3.set_row(0, points_0[i].x() * M0.get_row(2) - M0.get_row(0));
            A3.set_row(1, points_0[i].y() * M0.get_row(2) - M0.get_row(1));
            A3.set_row(2, points_1[i].x() * M3.get_row(2) - M3.get_row(0));
            A3.set_row(3, points_1[i].y() * M3.get_row(2) - M3.get_row(1));

            A4.set_row(0, points_0[i].x() * M0.get_row(2) - M0.get_row(0));
            A4.set_row(1, points_0[i].y() * M0.get_row(2) - M0.get_row(1));
            A4.set_row(2, points_1[i].x() * M4.get_row(2) - M4.get_row(0));
            A4.set_row(3, points_1[i].y() * M4.get_row(2) - M4.get_row(1));

            // STEP 4. Compute the SVD of A.
            Matrix44 H1, I1, J1;
            svd_decompose(A1,H1,I1,J1);
            Vector4D P1 = J1.get_column(J1.cols() - 1);
            // The 3D point P is in homogenous coordinates: a 3D point is represented as a 4D vector where the last element is a scale factor.
            P1 /= P1[3]; //divide by last element of P to scale
            P1_3D = P1.cartesian();
            P1_array.push_back(P1_3D);
            //std:: cout << "P1: " << P1_3D << std::endl;

            Matrix44 H2, I2, J2;
            svd_decompose(A2,H2,I2,J2);
            Vector4D P2 = J2.get_column(J2.cols() - 1);
            P2 /= P2[3]; //divide by last element of P to scale
            P2_3D = P2.cartesian();
            P2_array.push_back(P2_3D);
            //std:: cout << "P2: " << P2_3D << std::endl;

            Matrix44 H3, I3, J3;
            svd_decompose(A3,H3,I3,J3);
            Vector4D P3 = J3.get_column(J3.cols() - 1);
            P3 /= P3[3]; //divide by last element of P to scale
            P3_3D = P3.cartesian();
            P3_array.push_back(P3_3D);
            //std:: cout << "P3: " << P3_3D << std::endl;

            Matrix44 H4, I4, J4;
            svd_decompose(A4,H4,I4,J4);
            Vector4D P4 = J4.get_column(J4.cols() - 1);
            P4 /= P4[3]; //divide by last element of P to scale
            P4_3D = P4.cartesian();
            P4_array.push_back(P4_3D);
            //std:: cout << "P4: " << P4_3D << std::endl;
        }

        //STEP 5. Write recovered 3D points into 'points_3d':
        // Select the P with more positive zs value

        int zs_P1 = 0;

        for (int i = 0; i < points_0.size(); i++) {
            if (P1_array[i].z() > 0 && (Rt1*P1_array[i].homogeneous())[2]>0) {
                zs_P1 = zs_P1 + 1;
            } else if (P1_array[i].z() < 0) {
                zs_P1 = zs_P1;
            }
        }
        std::cout << "zs_P1: " << zs_P1 << std::endl;

        int zs_P2 = 0;

        for (int i = 0; i < points_0.size(); i++) {
            if (P2_array[i].z() > 0  && (Rt2*P2_array[i].homogeneous())[2]>0) {
                zs_P2 = zs_P2 + 1;
            } else if (P2_array[i].z() < 0) {
                zs_P2 = zs_P2;
            }
        }
        std::cout << "zs_P2: " << zs_P2 << std::endl;

        int zs_P3 = 0;

        for (int i = 0; i < points_0.size(); i++) {
            if (P3_array[i].z() > 0 && (Rt3*P3_array[i].homogeneous())[2]>0) {
                zs_P3 = zs_P3 + 1;
            } else if (P3_array[i].z() < 0) {
                zs_P3 = zs_P3;
            }
        }
        std::cout << "zs_P3: " << zs_P3 << std::endl;

        int zs_P4 = 0;

        for (int i = 0; i < points_0.size(); i++) {
            if (P4_array[i].z() > 0 && (Rt4*P4_array[i].homogeneous())[2]>0) {
                zs_P4 = zs_P4 + 1;
            } else if (P4_array[i].z() < 0) {
                zs_P4 = zs_P4;
            }
        }
        std::cout << "zs_P4: " << zs_P4 << std::endl;

        /*std::vector<Vector3D> P_array;
        if (zs_P1  > zs_P2 && zs_P1 > zs_P3 && zs_P1 > zs_P4) {
            std::cout << "R1 and t1 are the optimal" << std::endl;
            P_array = P1_array;
        }
        else if (zs_P2  > zs_P1 && zs_P2 > zs_P3 && zs_P1 > zs_P4) {
            std::cout << "R1 and t2 are the optimal" << std::endl;
            P_array = P2_array;
        }
        else if (zs_P3  > zs_P1 && zs_P3 > zs_P2 && zs_P3 > zs_P4) {
            std::cout << "R2 and t1 are the optimal" << std::endl;
            P_array = P3_array;
        }
        else if (zs_P4  > zs_P1 && zs_P4 > zs_P2 && zs_P4 > zs_P3) {
            std::cout << "R2 and t2 are the optimal" << std::endl;
            P_array = P4_array;
        }*/


        // Need to populate points_3d with variable which represents a collection of Vector3D objects.
        points_3d = P2_array;
        R = R1;
        t = t2;


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
        std::cout<<points_3d.size();
        for (size_t i = 0; i < points_3d.size(); i++)
        {
            std::cout<<points_3d[i]<<std::endl;
        }

        return points_3d.size() > 0;

    }
}