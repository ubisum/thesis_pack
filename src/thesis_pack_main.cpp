#include "stdio.h"
#include "g2o/core/factory.h"
#include "misc_utilities.h"
#include "matching_utilities.h"
#include "line_extractor.h"
#include "least_squares_utilities.h"
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
            float first = coord.first;
            float second = coord.second;
            temp_vec << first, second;
            //temp_vec << coord.first, coord.second;

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

    printLines(extractedLines[0], 0);
    printLines(extractedLines[1], 1);
    printLines(extractedLines[2], 2);


    MatrixXf prova(8,10);
    Vector4f vec1(4,1);
    vec1 << 1,2,3,4;
    Vector4f vec2(4,1);
    vec2 << 8,9,10,11;

    prova.block(0,0,4,10) = vec1.replicate(1,10);
    prova.block(4,0,4,10) = vec2.replicate(1,10);

    cout << lsIteration(T,prova,0.5,0.5) << endl;
    // cout << computeError(0,T, prova);



}
