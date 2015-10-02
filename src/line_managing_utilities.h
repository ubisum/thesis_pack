#ifndef LINE_MANAGING_UTILITIES_H
#define LINE_MANAGING_UTILITIES_H

#include "least_squares_utilities.h"
#include "misc_utilities.h"
#define NORM_THRESHOLD 0.1*0.1
#define MP_THRESHOLD 0.1*0.1

using namespace least_squares;
using namespace utilities;

namespace utilities
{

Vector2f middlepoint(Vector2f point_a, Vector2f point_b)
{
    Vector2f mp;
    mp << point_a(0)-point_b(0), point_a(1)-point_b(1);

    return 0.5*mp;
}

Vector2f projectPoint(Vector4f line_extremes, Vector2f point)
{
    // isolate coordinates
    float x1 = line_extremes(0);
    float y1 = line_extremes(1);
    float x2 = line_extremes(2);
    float y2 = line_extremes(3);

    // compute coefficients m and b
    float m = -(y2-y1)/(x1-x2);
    float b = y1-m*x1;

    // projection
    Vector2f projection;
    projection << (m*point(1) + point(0) - m*b) / (m*m + 1),
                  (m*m*point(1) + m*point(0) + b) / (m*m + 1);

    // return projection
    return projection;
}

int getPointPosition(Vector4f line_extremes, Vector2f projected_point)
{
    // isolate coordinates
    float xa = line_extremes(0);
    float xb = line_extremes(2);

    // compute parameter t
    float t = (projected_point(0) - xa)/(xb-xa);

    // return right value
    if (t<0)
        return -1;
    else if (t>=0 && t<=1)
        return 0;
    else
        return 1;
}

MatrixXf obtainNewExtremes(Vector4f line_left, Vector4f line_right)
{
    // isolate right extremes
    Vector2f right_ex1 = line_right.block(0,0,2,1);
    Vector2f right_ex2 = line_right.block(2,0,2,1);

    // positions of extremes
    int pos_ex1 = getPointPosition(line_left,right_ex1);
    int pos_ex2 = getPointPosition(line_left, right_ex2);

    // output matrix
    MatrixXf output;

    if((pos_ex1 == 1 && pos_ex2 == 1) || (pos_ex1 == -1 && pos_ex2 == -1))
    {
        output = MatrixXf::Zero(4,2);
        output.block(0,0,4,1) = line_left;
        output.block(0,1,4,1) = line_right;
    }

    else
    {
        output = MatrixXf::Zero(4,1);
        output.block(0,0,4,1) = line_left;

        if(pos_ex1 == -1)
            output.block(0,0,2,1) = right_ex1;

        else if(pos_ex1 == 1)
            output.block(2,0,2,1) = right_ex1;

        if(pos_ex2 == -1)
            output.block(0,0,2,1) = right_ex2;

        else if(pos_ex2 == 1)
            output.block(2,0,2,1) = right_ex2;
    }

    return output;
}

MatrixXf removeColumn(MatrixXf mat, int index)
{
    // output matrix
    MatrixXf output = mat;

    if(index < 0 || index > mat.cols()-1)
    {
        cout << "Error: the column with index " << index << " does not exist. The matrix has " <<
                mat.cols() << " columns." << endl << endl;
        return output;
    }

    if(mat.cols() > 1)
    {
        output = MatrixXf::Zero(mat.rows(), mat.cols()-1);
        cout << output << endl;

        for(int i = 0; i<index; i++)
            output.block(0,i,mat.rows(),1) = mat.block(0,i,mat.rows(),1);

        for(int j = index; j<output.cols(); j++)
            output.block(0,j,mat.rows(),1) = mat.block(0,j+1,mat.rows(),1);
    }

    return output;
}

MatrixXf mergeLines(MatrixXf left_lines, MatrixXf right_lines, MatrixXf transf_mat)
{
    // constants
    MatrixXf copy_left = left_lines;
    MatrixXf output;

    for(int i = 0; i<right_lines.cols(); i++)
    {
        // right normal
        Vector2f right_normal= right_lines.block(2,i,2,1);

        // transform right normal
        Vector2f transf_right_normal = transf_mat.block(0,0,2,2)*right_normal;

        // right middlepoint
        Vector2f right_middlepoint = right_lines.block(0,i,2,1);

        // transform right middlepoint
        Vector2f rx_mp_transf = transformPoints(right_middlepoint, transf_mat);

        // transform right extremes
        Vector4f transf_rx = Vector4f::Zero(4,1);
        transf_rx.block(0,0,2,1) = transformPoints(right_lines.block(6,i,2,1), transf_mat);
        transf_rx.block(2,0,2,1) = transformPoints(right_lines.block(8,i,2,1), transf_mat);

        // merged lines
        MatrixXf merged_lines;
        int column_index = -1;

        for(int j = 0; j<copy_left.cols(); j++)
        {
            // left normal
            Vector2f left_normal = copy_left.block(2,j,2,1);

            // difference between normals
            Vector2f normal_diff = left_normal-transf_right_normal;

            // difference between middlepoints
            Vector2f mp_diff = rx_mp_transf - copy_left.block(0,j,2,1);

            if(normal_diff.squaredNorm() < NORM_THRESHOLD && mp_diff.squaredNorm() < MP_THRESHOLD)
            {
                // new extremes
                merged_lines = obtainNewExtremes(copy_left.block(6,j,4,1), transf_rx);

                // save index
                column_index = j;

                // exit loop
                break;
            }
        }

        // add merged lines
        if(merged_lines.cols() > 0)
        {
            // new lines
            MatrixXf matrix_to_add;

            if(merged_lines.cols() == 1)
            {
                // create a new column
                MatrixXf new_line = MatrixXf::Zero(10,1);

                // fill column
                new_line.block(0,0,2,1) = middlepoint(merged_lines.block(0,0,2,1), merged_lines.block(2,0,2,1));
                new_line.block(2,0,4,1) = copy_left.block(2,column_index,4,1);
                new_line.block(4,0,2,1) = left_lines.block(4,column_index,2,1);
                new_line.block(6,0,4,1) = merged_lines;

                // save
                matrix_to_add = new_line;
            }

            else if(merged_lines.cols() == 2)
            {
                // create a new column
                MatrixXf new_line = MatrixXf::Zero(10,1);

                // fill column
                new_line.block(0,0,2,1) = rx_mp_transf;
                new_line.block(2,0,2,1) = transf_right_normal;
                new_line.block(4,0,2,1) = transformRT(right_lines.block(4,i,2,1), transf_mat);
                new_line.block(6,0,4,1) = transf_rx;

                // update
                merged_lines.block(0,1,10,1) = new_line;

                // save
                matrix_to_add = merged_lines;
            }

            // add lines
            MatrixXf updated_lines = MatrixXf::Zero(output.rows(), output.cols() + matrix_to_add.cols());
            updated_lines.block(0,0,output.rows(), output.cols()) = output;
            updated_lines.block(0, output.cols(), output.rows(), matrix_to_add.cols()) = matrix_to_add;
            output = updated_lines;

            // remove column
            removeColumn(copy_left,column_index);

        }

        // no match in left scan
        else
        {
            // create a new column
            MatrixXf new_line = MatrixXf::Zero(10,1);

            // fill column
            new_line.block(0,0,2,1) = rx_mp_transf;
            new_line.block(2,0,2,1) = transf_right_normal;
            new_line.block(4,0,2,1) = transformRT(right_lines.block(4,i,2,1), transf_mat);
            new_line.block(6,0,4,1) = transf_rx;

            // create a new matrix of lines
            MatrixXf updated_lines = MatrixXf::Zero(output.rows(), output.cols() + 1);
            updated_lines.block(0,0,output.rows(), output.cols()) = output;
            updated_lines.block(0,output.cols(), 10, 1) = new_line;

            // save
            output = updated_lines;

        }

    }

    // return lines
    return output;
}

}

#endif // LINE_MANAGING_UTILITIES_H
