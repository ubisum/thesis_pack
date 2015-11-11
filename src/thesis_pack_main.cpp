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

    Vector4f a(15,12,24,20);
    obtainNewExtremes(a,Vector4f(15,16,19,22));
    cout << "Projection " << endl << projectPoint(Vector4f(15,12,24,20), Vector2f(10,14)) << endl;

//    vector<vecPairsList> reducedVector(extractedLines.begin()+1, extractedLines.begin()+2);
    int init = 0;
    int end = 1;
    vector<vecPairsList> reducedVector;
    for(int k = init; k<= end; k++)
        reducedVector.push_back(extractedLines[k]);
    //MatrixXf m = mergeLines(reducedVector);
    MatrixXf m = mergeLines(extractedLines);
    cout << "Size m: " << m.rows() << " " << m.cols() << endl;
    printLinesByExtremes(m.block(6,0,4,m.cols()), MatrixXf::Identity(3,3), "merged_lines.txt");

//    cout << "prove" << endl << endl;
//    cout << obtainNewExtremes(Vector4f(7,7,10,10), Vector4f(3,3,5,5)) << endl << endl;
//    cout << obtainNewExtremes(Vector4f(7,7,10,10), Vector4f(12,12,14,14)) << endl << endl;
//    cout << obtainNewExtremes(Vector4f(7,7,10,10), Vector4f(5,5,8,8)) << endl << endl;
//    cout << obtainNewExtremes(Vector4f(7,7,10,10), Vector4f(8,8,11,11)) << endl << endl;

}
