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

vector<vecPairsList> createVector()
{
    vecPair left1(Vector2f(2,1), Vector2f(5,1));
    vecPair left2(Vector2f(2,1), Vector2f(2,3));

//    vecPair right1(Vector2f(2.01,1.01), Vector2f(5.01,1.01));
//    vecPair right2(Vector2f(2.01,1.01), Vector2f(2.01,3.01));

        vecPair right1(Vector2f(3,2), Vector2f(6,2));
        vecPair right2(Vector2f(3,2), Vector2f(3,4));

    vecPairsList vpl1;
    vpl1.push_back(left1);
    vpl1.push_back(left2);

    vecPairsList vpl2;
    vpl2.push_back(right1);
    vpl2.push_back(right2);

    //return vector<vecPairsList>(vpl1, vpl2);
    vector<vecPairsList> output;
    output.push_back(vpl1);
    output.push_back(vpl2);

    return output;


}

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
    //cout << "Projection " << endl << projectPoint(Vector4f(15,12,24,20), Vector2f(10,14)) << endl;

     cout << "Projection " << endl << projectByEquation(Vector3f(2,-1,0), Vector2f(2,3)) << endl;

//    vector<vecPairsList> reducedVector(extractedLines.begin()+1, extractedLines.begin()+2);
    int init = 0;
    int end = 1;
    vector<vecPairsList> reducedVector;
    for(int k = init; k<= end; k++)
        reducedVector.push_back(extractedLines[k]);
    //MatrixXf m = mergeLines(createVector());
    MatrixXf m = mergeLines(reducedVector);
    //MatrixXf m = mergeLines(extractedLines);
    cout << "Size m: " << m.rows() << " " << m.cols() << endl;
    printLinesByExtremes(m.block(6,0,4,m.cols()), MatrixXf::Identity(3,3), "merged_lines.txt");

    cout << endl << "Proiezione punto" << endl;
    cout << projectByEquation(Vector3f(2, -1, 0), Vector2f(2,3)) << endl << endl;


    Vector4f line1(2,3,2,1);
    Vector4f line2(2,1,5,1);
    MatrixXf two_lines(4,2);
    two_lines.block(0,0,4,1) = line1;
    two_lines.block(0,1,4,1) = line2;

    float angle = 45*M_PI/180;
    MatrixXf prova_t(3,3);
    prova_t << cos(angle), -sin(angle), 5, sin(angle), cos(angle), 2, 0, 0, 1;
//    prova_t << 1,0,5,0,1,0,0,0,1;
//    MatrixXf prova_t = MatrixXf::Identity(3,3);

    MatrixXf transf_vectors(4,2);
    transf_vectors.block(0,0,2,2) = transformVectors(two_lines.block(0,0,2,2), prova_t);
    transf_vectors.block(2,0,2,2) = transformVectors(two_lines.block(2,0,2,2), prova_t);
    cout << "transf_vectors" << endl << transf_vectors << endl;

    remove("prova_transf.txt");
    FILE* prova_transf = fopen("prova_transf.txt", "a");
    stringstream ss_prova;
    for(int w = 0; w<transf_vectors.cols(); w++)
    {
        Vector4f column = transf_vectors.block(0,w,4,1);

        ss_prova << column(0) << "\t" << column(1) << "\n" <<
                    column(2) << "\t" << column(3) << "\n\n";

        fputs(ss_prova.str().c_str(), prova_transf);
    }

    Vector2f polar1 = polarRepresentation(line1.block(0,0,2,1), line1.block(2,0,2,1));
    //Vector2f polar2 = polarRepresentation(transf_vectors.block(0,0,2,1), transf_vectors.block(2,0,2,1));
    Vector2f polar2 = new_transformRT(transf_vectors.block(0,0,4,1));

    cout << "differenza theta " << 180*fabs(polar1(1)-polar2(1))/M_PI << endl;


}
