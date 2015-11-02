#ifndef LINE_MANAGING_UTILITIES_H
#define LINE_MANAGING_UTILITIES_H

#include "least_squares_utilities.h"
#include "misc_utilities.h"

#define NORM_THRESHOLD 0.8
#define MP_THRESHOLD pow(10,-3)
#define ANGLE_THRESHOLD 5*M_PI/180
#define RHO_THRESHOLD 0.1

using namespace least_squares;
using namespace utilities;

namespace utilities
{

typedef pair<int, float> norm_index;
typedef pair<int, int> indeces_pair;

Vector2f middlepoint(Vector2f point_a, Vector2f point_b)
{
    Vector2f mp;
    mp << point_a(0)+point_b(0), point_a(1)+point_b(1);

    return 0.5*mp;
}

norm_index findClosestNorm(MatrixXf matrix_of_normals, Vector2f norm)
{
    if(matrix_of_normals.cols() == 0)
        return norm_index(-1,-1);

    else
    {
        int min_index = 0;
        float min_sq_normal = (matrix_of_normals.block(0,0,2,1)-norm).squaredNorm();

        for(int i = 1; i<matrix_of_normals.cols(); i++)
        {
            float temp_diff = (matrix_of_normals.block(0,i,2,1)-norm).squaredNorm();
            if(temp_diff < min_sq_normal)
            {
                min_index = i;
                min_sq_normal = temp_diff;
            }
        }

        return norm_index(min_index, min_sq_normal);
    }

}

Vector2f projectPoint(Vector4f line_extremes, Vector2f point)
{
    // isolate coordinates
    float x1 = line_extremes(0);
    float y1 = line_extremes(1);
    float x2 = line_extremes(2);
    float y2 = line_extremes(3);

    Vector2f ab = -(Vector2f(x1, y1)-Vector2f(x2, y2));
    Vector2f projection = point.dot(ab)/ ab.dot(ab) * ab;

    // return projection
    return projection;
}

int getPointPosition(Vector4f line_extremes, Vector2f point)
{
    // project point
    Vector2f projected_point = projectPoint(line_extremes,point);

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
    int pos_ex1 = getPointPosition(line_left, right_ex1);
    int pos_ex2 = getPointPosition(line_left, right_ex2);

    // output matrix
    MatrixXf output;

    if((pos_ex1 == 1 && pos_ex2 == 1) || (pos_ex1 == -1 && pos_ex2 == -1))
    {
        //cout << "Caso 1" << endl;
        output = MatrixXf::Zero(4,2);
        output.block(0,0,4,1) = line_left;
        output.block(0,1,4,1) = line_right;
    }

    else if(pos_ex1 == 0 && pos_ex2 == 0)
    {
        //cout << "Caso 2" << endl;
        output = MatrixXf::Zero(4,1);
        output = line_left;
    }

    else
    {
        //cout << "Merging lines" << endl;
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

        for(int i = 0; i<index; i++)
            output.block(0,i,mat.rows(),1) = mat.block(0,i,mat.rows(),1);

        for(int j = index; j<output.cols(); j++)
            output.block(0,j,mat.rows(),1) = mat.block(0,j+1,mat.rows(),1);
    }

    return output;
}

vector<indeces_pair> mergingIteration(MatrixXf left_lines, MatrixXf right_lines, MatrixXf T)
{
    remove("mergeIteration.txt");
    remove("mp_lines.txt");
    FILE* merge_file = fopen("mergeIteration.txt", "a");
    FILE* mp_file = fopen("mp_lines.txt", "a");
    vector<indeces_pair> pairs_list;

    for(int i = 0; i<right_lines.cols(); i++)
    {
        // rotate normal
        Vector2f rot_normal = T.block(0,0,2,2)*right_lines.block(2,i,2,1);

        for(int j = 0; j<left_lines.cols(); j++)
        {
            // left normal
            Vector2f left_normal = left_lines.block(2,j,2,1);

            if(rot_normal.transpose()*left_normal > NORM_THRESHOLD)
            {
                stringstream ss_lines, ss_mp;

                // transform points
                Vector2f mp_rot = transformVectors(right_lines.block(0,i,2,1), T);
                Vector2f ex1_rot = transformVectors(right_lines.block(6,i,2,1), T);
                Vector2f ex2_rot = transformVectors(right_lines.block(8,i,2,1), T);

                // write middlepoints to file
                Vector2f left_mp = left_lines.block(0,j,2,1);
                ss_mp << left_mp(0) << "\t" << left_mp(1) << "\n" <<
                           mp_rot(0) << "\t" << mp_rot(1) << "\n\n";
                fputs(ss_mp.str().c_str(), mp_file);

                // write lines to file
                Vector2f left_ex1 = left_lines.block(6,j,2,1);
                Vector2f left_ex2 = left_lines.block(8,j,2,1);

                ss_lines << left_ex1(0) << "\t" << left_ex1(1) << "\n" <<
                            left_ex2(0) << "\t" << left_ex2(1) << "\n\n" <<
                            ex1_rot(0) << "\t" << ex1_rot(1) << "\n" <<
                            ex2_rot(0) << "\t" << ex2_rot(1) << "\n\n";
                fputs(ss_lines.str().c_str(), merge_file);

                pairs_list.push_back(indeces_pair(j, i));
            }
        }
    }

    fclose(merge_file);
    fclose(mp_file);

    return pairs_list;
}

void evaluateByRho(vector<indeces_pair> indeces, MatrixXf left_lines, MatrixXf right_lines, MatrixXf T)
{
    remove("evaluate_rho.txt");
    FILE* rho_fid = fopen("evaluate_rho.txt", "a");

    for(int i = 0; i<(int)indeces.size(); i++)
    {
        // get a pair
        indeces_pair ip = indeces[i];
        cout << ip.first << " " << ip.second << endl;

        // get left rho
        float left_rho = left_lines(4, ip.first);

        // get right rho
        Vector2f rt_rot = transformRT(right_lines.block(4,ip.second,2,1), T);
        float right_rho = rt_rot(0);

        if(fabs(right_rho-left_rho) < RHO_THRESHOLD)
        {
            stringstream ss_rho;

            // transform right middlepoint
            Vector2f mp_transf = transformVectors(right_lines.block(0, ip.second, 2, 1), T);
            Vector2f left_mp = left_lines.block(0,ip.first,2,1);
            // write to file
            ss_rho << left_mp(0) << "\t" << left_mp(1) << "\n" <<
                      mp_transf(0) << "\t" << mp_transf(1) << "\n\n";

            fputs(ss_rho.str().c_str(), rho_fid);
        }
    }

    fclose(rho_fid);
}

MatrixXf mergeLines(const vector<vecPairsList>& extractedLines)
{

    MatrixXf final_lines;

    if(extractedLines.size() == 1)
        final_lines = linesByCol(extractedLines,0);

    else if (extractedLines.size() > 1)
    {
        // final lines
        final_lines = linesByCol(extractedLines,0);
        printLinesByExtremes(final_lines.block(6,0,4, final_lines.cols()), MatrixXf::Identity(3,3), "convertedLines_0.txt");
        for(int i = 1; i<(int)extractedLines.size(); i++)
        {
//            if(i < 2)
//            {
                // get lines in matrix form
                MatrixXf current_lines = linesByCol(extractedLines,i);

                // get current transformation matrix
                MatrixXf out_li, out_lj;
                MatrixXf T = getTransformationMatrix(final_lines, current_lines, out_li, out_lj);
                printLinesByExtremes(current_lines.block(6,0,4,current_lines.cols()), T, "convertedLines_1.txt");

                // merge current lines
                vector<indeces_pair> pairs = mergingIteration(final_lines, current_lines, T);
                evaluateByRho(pairs, final_lines, current_lines, T);
                //cout << merged_lines.rows() << " " << merged_lines.cols() << endl;
                //cout << final_lines.rows() << " " << final_lines.cols() << endl;
                cout << "Dopo Merged Lines" << endl;

                // update
//                final_lines = MatrixXf::Zero(merged_lines.rows(), merged_lines.cols());
//                final_lines.block(0,0,merged_lines.rows(), merged_lines.cols()) = merged_lines;
                //final_lines = merged_lines;
//            }

        }
    }

    // return final lines
    return final_lines;
}

}

#endif // LINE_MANAGING_UTILITIES_H
