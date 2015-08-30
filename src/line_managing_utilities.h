#ifndef LINE_MANAGING_UTILITIES_H
#define LINE_MANAGING_UTILITIES_H

#include "least_squares_utilities.h"
#include "misc_utilities.h"

using namespace least_squares;
using namespace utilities;

namespace utilities
{

MatrixXf mergeLines(MatrixXf left_lines, MatrixXf right_lines, MatrixXf Li, MatrixXf Lj, MatrixXf T)
{
    // create larger matrices
    MatrixXf larger_left = MatrixXf::Zero(left_lines.rows()+1, left_lines.cols());
    MatrixXf larger_right = MatrixXf::Zero(right_lines.rows()+1, right_lines.cols());

    // fill larger matrices
    larger_left.block(0,0,left_lines.rows(),left_lines.cols()) = left_lines;
    larger_right.block(0,0,right_lines.rows(), right_lines.cols()) = right_lines;

    // insert flags
    for(int i = 0; i<Li.cols(); i++)
        larger_left(left_lines.rows(),Li(10,i)) = 1;

    for(int j = 0; j<Lj.cols(); j++)
        larger_right(right_lines.rows(), Lj(10,j)) = 1;

    // count valid lines
    MatrixXf valid_left = larger_left.block(left_lines.rows(), 0, 1, left_lines.cols())*
            MatrixXf::Constant(left_lines.cols(),1,1);

    MatrixXf valid_right = larger_right.block(right_lines.rows(), 0, 1, right_lines.cols())*
            MatrixXf::Constant(right_lines.cols(),1,1);

    // create a new matrix
    MatrixXf valid_lines = MatrixXf::Zero(10, valid_left(0,0) + right_lines.cols()-valid_right(0,0));

    // fill new matrix
    int counter = 0;

    for(int counterLeft = 0; counterLeft<larger_left.cols(); counterLeft++)
    {
        if(larger_left(10,counterLeft) == 1)
        {
            valid_lines.block(0,counter,10,1) = larger_left.block(0,counterLeft,10,1);
            counter++;
        }
    }

    for(int counterRight = 0; counterRight<larger_right.cols(); counterRight++)
    {
        if(larger_right(10,counterRight) == 0)
        {
            // get a right column
            MatrixXf temp_col = larger_right.block(0, counterRight,10,1);

            // fill current new column
            valid_lines.block(0,counter,4,1) = transformPoints(T, temp_col.block(0,0,4,1));
            valid_lines.block(4,counter,2,1) = transformRT(temp_col.block(4,0,2,1), T);
            valid_lines.block(6,counter,2,1) = transformVectors(temp_col.block(6,0,2,1), T);
            valid_lines.block(8,counter,2,1) = transformVectors(temp_col.block(8,0,2,1), T);

            // increment counter
            counter++;
        }
    }

    // return valid lines
    return valid_lines;

}

}

#endif // LINE_MANAGING_UTILITIES_H
