#ifndef LINE_MANAGING_UTILITIES_H
#define LINE_MANAGING_UTILITIES_H

#include "least_squares_utilities.h"
#include "misc_utilities.h"

#define NORM_THRESHOLD 0.8
#define MP_THRESHOLD pow(0.3,2)
#define ANGLE_THRESHOLD 5*M_PI/180
#define THETA_THRESHOLD 0.1

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

Vector2f middlepoint(Vector4f line_extremes)
{
    return middlepoint(line_extremes.block(0,0,2,1), line_extremes.block(2,0,2,1));
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

void addRemainingCols(MatrixXf& final_lines, MatrixXf right_lines, vector<indeces_pair> pairs, MatrixXf T)
{
    // create enlarged matrix
    MatrixXf enlarged = MatrixXf::Zero(11, right_lines.cols());
    enlarged.block(0,0,right_lines.rows(),right_lines.cols()) = right_lines;

    // set booleans
    for(int i = 0; i<(int)pairs.size(); i++)
        enlarged(10, pairs[i].second) = 1;

    // add columns
    for(int j = 0; j<enlarged.cols(); j++)
    {
        if(enlarged(10,j) == 0)
        {
            // get and transform column
            MatrixXf column = transformFullLine(enlarged.block(0, j, 10, 1), T);

            // create a new matrix
            MatrixXf new_Matrix = MatrixXf::Zero(final_lines.rows(), final_lines.cols()+1);

            // copy data
            new_Matrix.block(0,0,final_lines.rows(), final_lines.cols()) = final_lines;

            // add columns
            new_Matrix.block(0,final_lines.cols(), 10, 1) = column;

            // update
            final_lines = new_Matrix;
        }
    }


}

Vector2f projectPoint(Vector4f line_extremes, Vector2f point)
{
    // isolate coordinates
    float x1 = line_extremes(0);
    float y1 = line_extremes(1);
    float x2 = line_extremes(2);
    float y2 = line_extremes(3);

//    Vector2f ab = -(Vector2f(x1, y1)-Vector2f(x2, y2));
//    Vector2f projection = point.dot(ab)/ ab.dot(ab) * ab;

//    // return projection
//    return projection;

    MatrixXf A(2,2);
    A << x2-x1, y2-y1, y1-y2, x2-x1;

    Vector2f b;
    b << -point(0)*(x2-x1)-point(1)*(y2-y1), -y1*(x2-x1)+x1*(y2-y1);

    return A.colPivHouseholderQr().solve(-b);
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

    // projections of points
    Vector2f proj_1 = projectPoint(line_left, right_ex1);
    Vector2f proj_2 = projectPoint(line_left, right_ex2);

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
//            output.block(0,0,2,1) = right_ex1;
            output.block(0,0,2,1) = proj_1;

        else if(pos_ex1 == 1)
//            output.block(2,0,2,1) = right_ex1;
            output.block(2,0,2,1) = proj_1;


        if(pos_ex2 == -1)
//            output.block(0,0,2,1) = right_ex2;
            output.block(0,0,2,1) = proj_2;


        else if(pos_ex2 == 1)
//            output.block(2,0,2,1) = right_ex2;
            output.block(2,0,2,1) = proj_2;

    }

    remove("extremes.txt");
    remove("projections.txt");

    FILE* extremes = fopen("extremes.txt", "a");
    FILE* projections = fopen("projections.txt", "a");

    stringstream ss_extremes;
    ss_extremes << line_left(0) << "\t" << line_left(1) << "\n" <<
                   line_left(2) << "\t" << line_left(3);
    fputs(ss_extremes.str().c_str(), extremes);

    for(int i = 0; i<output.cols(); i++)
    {
        stringstream ss_projections;
        ss_projections << output(0,i) << "\t" << output(1,i) << "\n" <<
                          output(2,i) << "\t" << output(3,i) << "\n\n";

        fputs(ss_projections.str().c_str(), projections);
    }

    fclose(extremes);
    fclose(projections);

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

vector<indeces_pair> evaluateByNormal(MatrixXf left_lines, MatrixXf right_lines, MatrixXf T)
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

            // middlepoints
            Vector2f left_mp = left_lines.block(0,0,2,1);
            Vector2f mp_rot = transformVectors(right_lines.block(0,i,2,1), T);

            if(/*rot_normal.transpose()*left_normal*/ rot_normal.dot(left_normal) > NORM_THRESHOLD
                   /*&& left_mp.transpose()*mp_rot < MP_THRESHOLD*/)
            {
                stringstream ss_lines, ss_mp;

                // transform points
                //Vector2f mp_rot = transformVectors(right_lines.block(0,i,2,1), T);
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

vector<indeces_pair> evaluateByTheta(vector<indeces_pair> indeces, MatrixXf left_lines, MatrixXf right_lines, MatrixXf T)
{
    remove("evaluate_rho.txt");
    FILE* rho_fid = fopen("evaluate_rho.txt", "a");
    vector<indeces_pair> after_evaluation;

    for(int i = 0; i<(int)indeces.size(); i++)
    {
        // get a pair
        indeces_pair ip = indeces[i];

        // get left theta
        float left_theta = left_lines(5, ip.first);

        // middlepoints
        Vector2f left_mp = left_lines.block(0,ip.first,2,1);
        Vector2f mp_rot = transformVectors(right_lines.block(0,ip.second,2,1), T);

        // get right theta
        Vector2f rt_rot = transformRT(right_lines.block(4,ip.second,2,1), T);
        float right_theta = rt_rot(1);

        if(fabs(right_theta-left_theta) < THETA_THRESHOLD && (left_mp-mp_rot).transpose()*(left_mp-mp_rot) < MP_THRESHOLD)
        {
            stringstream ss_rho;

            // transform right middlepoint
            Vector2f mp_transf = transformVectors(right_lines.block(0, ip.second, 2, 1), T);
            Vector2f left_mp = left_lines.block(0,ip.first,2,1);
            // write to file
            ss_rho << left_mp(0) << "\t" << left_mp(1) << "\n" <<
                      mp_transf(0) << "\t" << mp_transf(1) << "\n\n";

            fputs(ss_rho.str().c_str(), rho_fid);
            after_evaluation.push_back(ip);
        }
    }

    fclose(rho_fid);
    return after_evaluation;
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
            // get lines in matrix form
            MatrixXf current_lines = linesByCol(extractedLines,i);
            cout << "Elaborating scan " << i+1 << "/" << extractedLines.size() << endl;

            // get current transformation matrix
            MatrixXf out_li, out_lj;
            MatrixXf T = getTransformationMatrix(final_lines, current_lines, out_li, out_lj);
            printLinesByExtremes(current_lines.block(6,0,4,current_lines.cols()), T, "convertedLines_1.txt");

            // get merging pairs
            vector<indeces_pair> pairs = evaluateByNormal(final_lines, current_lines, T);
            cout << "Num normal: " << pairs.size() << endl;
            vector<indeces_pair> after_rho_evaluation = evaluateByTheta(pairs, final_lines, current_lines, T);
            cout << "Num eval: " << after_rho_evaluation.size() << " right lines: " << current_lines.cols() << endl;
            cout << "Lines to add: " << current_lines.cols()-(int)after_rho_evaluation.size() << endl;

            // loop over pairs
            for(int j = 0; j<(int)after_rho_evaluation.size(); j++)
            {
                // grab a pair
                indeces_pair ip = after_rho_evaluation[j];

                // grab left line
                Vector4f left_line = final_lines.block(6, ip.first, 4, 1);

                // grab right line
                Vector4f right_line = current_lines.block(6, ip.second, 4, 1);

                // transform right line
                Vector4f transf_right_line = Vector4f::Zero(4,1);
                transf_right_line.block(0,0,2,1) = transformVectors(right_line.block(0,0,2,1), T);
                transf_right_line.block(2,0,2,1) = transformVectors(right_line.block(2,0,2,1), T);

                // merge lines
                MatrixXf merged_lines = obtainNewExtremes(left_line, transf_right_line);

                // lines were merged
                if(merged_lines.cols() == 1)
                {
                    // create a new line column
                    MatrixXf new_col = MatrixXf::Zero(10,1);
                    new_col.block(0,0,2,1) = middlepoint(merged_lines);
                    new_col.block(2,0,4,1) = final_lines.block(2,ip.first,4,1);
                    new_col.block(6,0,4,1) = merged_lines;

                    // update column
                    final_lines.block(0,ip.first, 10, 1) = new_col.block(0,0,10,1);
                }

                // lines were not merged
                else
                {
                    // create a new line column
                    MatrixXf new_col = transformFullLine(current_lines.block(0,ip.second, 10, 1), T);

                    // create a larger matrix
                    MatrixXf new_lines = MatrixXf::Zero(final_lines.rows(), final_lines.cols()+1);

                    // copy data to new matrix
                    new_lines.block(0,0,final_lines.rows(), final_lines.cols()) = final_lines;

                    // insert new column
                    new_lines.block(0, final_lines.cols(), 10, 1) = new_col;

                    // update
                    final_lines = new_lines;
                }

            }

            // add other lines
            cout << "Prima add: " << final_lines.cols() << endl;
            addRemainingCols(final_lines, current_lines, after_rho_evaluation, T);
            cout << "Dopo add:" << final_lines.cols() << endl << endl;

        }
    }

    // return final lines
    //cout << "Size: " << final_lines.rows() << " " << final_lines.cols() << endl;
    return final_lines;
}

}

#endif // LINE_MANAGING_UTILITIES_H
