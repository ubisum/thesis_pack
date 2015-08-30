#include "stdio.h"
#include "g2o/core/factory.h"
#include "misc_utilities.h"
#include "matching_utilities.h"
#include "line_extractor.h"
#include "least_squares_utilities.h"
#include "line_managing_utilities.h"
#include "vector"

using namespace utilities;
using namespace least_squares;

int main(int argc, char** argv)
{
    // data
    vector<completeInformation> compl_info, selected_info;
    vector<vecPairsList> extractedLines;

    // read robot's info
    parseRobotInfo("/home/ubisum/fuerte_workspace/thesis_pack/src/files/robotInfo.txt", compl_info);

    // select relevant poses
    selected_info = selectFrames(compl_info);

    // compute lines for each pose
    for(int i = 0; i<(int)selected_info.size(); i++)
    {
        Vector2fVector vector;

        for(int j = 0; j<(int)selected_info[i].points.size(); j++)
        {
            // extract a coordinate of a scan's point
            coordinate coord = selected_info[i].points[j];

            // fill a vector
            Vector2f temp_vec;
            temp_vec << coord.first, coord.second;

            // store vector
            vector.push_back(temp_vec);
        }

        // compute and store lines for current scan
        extractedLines.push_back(computeLines(vector));
    }

    MatrixXf T(3,3);
        T <<  0.93969, -0.34202, 0.1,
                0.34202, 0.93969, 0,
            0, 0, 1;

    printLinesFullDescription(extractedLines,0);
    printLinesFullDescription(extractedLines,1);

    MatrixXf linesCol_0 = linesByCol(extractedLines,0);
    MatrixXf linesCol_1 = linesByCol(extractedLines, 1);

    MatrixXf out_li, out_lj;
    MatrixXf transf_mat =  getTransformationMatrix(linesCol_0, linesCol_1, out_li, out_lj);
    MatrixXf new_lines = mergeLines(linesCol_0, linesCol_1, out_li, out_lj, transf_mat);
    printLinesByExtremes(linesCol_0.block(6,0,4,linesCol_0.cols()), MatrixXf::Identity(3,3),"linesCol_0.txt");
    printLinesByExtremes(linesCol_1.block(6,0,4,linesCol_1.cols()), transf_mat, "linesCol_1.txt");
    printLinesByExtremes(new_lines.block(6,0,4,new_lines.cols()), MatrixXf::Identity(3,3), "new_lines.txt");
    printAssociations(out_li, out_lj, transf_mat, "assoc.txt");


    //printLines(linesCol_0.block(6,0,4,linesCol_0.cols()), MatrixXf::Identity(3,3), "convertedLines_0.txt");
    //printLines(linesCol_1.block(6,0,4,linesCol_1.cols()), MatrixXf::Identity(3,3), "convertedLines_1.txt");

    MatrixXf ne_i(6,linesCol_0.cols());
    ne_i.block(0,0,2,linesCol_0.cols()) = linesCol_0.block(2,0,2,linesCol_0.cols());
    ne_i.block(2,0,4,linesCol_0.cols()) = linesCol_0.block(6,0,4,linesCol_0.cols());
    MatrixXf ne_j(6,linesCol_1.cols());
    ne_j.block(0,0,2,linesCol_1.cols()) = linesCol_1.block(2,0,2,linesCol_1.cols());
    ne_j.block(2,0,4,linesCol_1.cols()) = linesCol_1.block(6,0,4,linesCol_1.cols());

//    //cout << linesCol_1 << endl;
    //cout << computeDistanceNM(ne_i,ne_j,0.1,0.4,T) << endl;
    //cout << ne_j << endl;
    //cout << linesCol_0.transpose() << endl;

    MatrixXf uno = MatrixXf::Zero(3,3);
    uno(0,0) = -1;
    uno(1,1) = -1;
    uno(2,0) = -1;
    uno(2,2) = -1;
    // cout << findMinIndeces(uno) << endl;
    //out << uno.minCoeff() << endl;
    MatrixXf dist = computeDistanceNM(ne_i, ne_j, 0.1, 0.4, T);
    MatrixXi assoc = computeAssociations(dist,0.8);
    //cout << assoc + MatrixXi::Constant(assoc.rows(), assoc.cols(),1) << endl << endl;
    //cout << computeAssociations(dist,0.8) << endl;
    //cout << findValueInMatrix(uno, uno.minCoeff()) << endl;
    //cout << findRepetitions(assoc) << endl;
    //findRepetitions(assoc);
    //cout << getTransformationMatrix(linesCol_0, linesCol_1) << endl << endl;
    //MatrixXf new_row = MatrixXf::Zero(10,17);
    //MatrixXf new_dist = addIndecesRow(dist);
    //cout << addIndecesRow(dist) << endl << endl;
    //cout << dist << endl;
    //cout << "originali " << dist.rows() << " " << dist.cols() << endl << endl;
    //cout << "dopo " << new_dist.rows() << " " << new_dist.cols() << endl << endl;
    //cout << new_row.block(new_row.rows()-1,0,1,new_row.cols()) << endl << endl;
    MatrixXf prova_rt(2,5);
    prova_rt << 1,2,3,4,5,6,7,8,9,10;
    cout << prova_rt << endl << endl;
    //cout << transformRT(prova_rt, T) << endl;

}
