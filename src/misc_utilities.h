#ifndef MISC_UTILITIES_H
#define MISC_UTILITIES_H

#include "matching_utilities.h"
#include "fstream"
#include "string"


#define MIN_LINE_LENGTH 0.1*0.1
#define MIN_FRAME_DIST 0.1*0.1

using namespace std;

namespace utilities
{

MatrixXf transformRT(MatrixXf lines, MatrixXf T)
{
    // local data
    float x = T(0,2);
    float y = T(1,2);
    float theta = atan2(T(1,0),T(0,0));

    // prepare output
    MatrixXf transRT = MatrixXf::Zero(lines.rows(), lines.cols());

    // transform
    for(int i = 0; i<lines.cols(); i++)
    {
        // select a line
        Vector2f line_ith = lines.block(0,i,2,1);

        // current values
        float rho_ith = line_ith(0);
        float theta_ith = line_ith(1);

        // compute transformation
        float new_angle = theta_ith + theta;
        transRT(1,i) = new_angle;
        transRT(0,i) = rho_ith + cos(new_angle)*x + sin(new_angle)*y;
    }

    // return transformed data
    return transRT;

}

MatrixXf transformVectors(MatrixXf vectors, MatrixXf T)
{
    // output matrix
    MatrixXf output_vectors = MatrixXf::Zero(vectors.rows()+1, vectors.cols());
    output_vectors.block(0,0,vectors.rows(),vectors.cols()) = vectors;

    // set last row in output matrix
    output_vectors.block(vectors.rows(),0,1,vectors.cols()) = MatrixXf::Constant(1, vectors.cols(), 1);

    // product
    MatrixXf product = T*output_vectors;

    // return product
    return product.block(0,0,vectors.rows(),vectors.cols());
}

void parseRobotInfo(string file_name, vector<completeInformation>& ci)
{
    ifstream file(file_name.c_str());
    string str;
    char* pch;
    vector<string> parseInfo;
    int infoIndex = -1;

    while(getline(file, str))
    {
        char buffer[str.length()];
        str.copy(buffer, str.length(), 0);
        pch = strtok (buffer," \t\n\r");
        while(pch!=NULL)
        {
            parseInfo.push_back(string(pch));
            pch = strtok (NULL," \t\n\r");
        }

        if(strcmp(parseInfo[0].c_str(), "POSITION") == 0)
        {
            completeInformation temp_ci;
            temp_ci.position = make_pair((float)atof(parseInfo[1].c_str()),
                                         (float)atof(parseInfo[2].c_str()));
            temp_ci.angle = (float)atof(parseInfo[3].c_str());

            ci.push_back(temp_ci);
            infoIndex++;
        }

        else if(strcmp(parseInfo[0].c_str(), "POINT") == 0)
        {
            coordinate temp_coord;
            temp_coord.first = (float)atof(parseInfo[1].c_str());
            temp_coord.second = (float)atof(parseInfo[2].c_str());

            ci[infoIndex].points.push_back(temp_coord);
        }


        parseInfo.clear();
    }

}

vector<completeInformation> selectFrames(const vector<completeInformation>& ci)
{
    vector<completeInformation> selectedFrames;

    if(ci.size() >0)
    {
        int index = 0;
        for(int i = 1; i<(int)ci.size(); i++)
        {
            completeInformation currentFrame = ci[i];
            completeInformation prevFrame = ci[index];

            Vector2f diffs;
            diffs << currentFrame.position.first - prevFrame.position.first,
                     currentFrame.position.second - prevFrame.position.second;

            if(diffs.transpose()*diffs >= MIN_FRAME_DIST || currentFrame.angle != prevFrame.angle)
            {
                selectedFrames.push_back(currentFrame);
                index = i;
            }
        }
    }

    return selectedFrames;
}

vector<vecPairsList> selectLines(const vector<vecPairsList> lines)
{
    vector<vecPairsList> selectedLines;

    for(int i = 0; i<(int)lines.size(); i++)
    {
        vecPairsList vpl = lines[i];
        vecPairsList new_vpl;

        for (int j = 0; j<(int)vpl.size(); j++)
        {
            vecPair vPair = vpl[j];
            Vector2f ex_1 = vPair.first;
            Vector2f ex_2 = vPair.second;

            Vector2f diffs;
            diffs << ex_1(0)-ex_2(0), ex_1(1)-ex_2(1);

            if(diffs.transpose()*diffs >= MIN_LINE_LENGTH)
                new_vpl.push_back(vPair);
        }

        selectedLines.push_back(new_vpl);
        //cout << "Remaining lines: " << new_vpl.size() << "/" << vpl.size() << endl;
    }

    return selectedLines;
}

void printLinesFullDescription(const vector<vecPairsList>& lines_vector, int index)
{
    //cout << "entrata" << endl;
    stringstream ss_filename;
    ss_filename << "convertedLines_" << index << ".txt";
    vecPairsList Lines = lines_vector[index];

    remove(ss_filename.str().c_str());
    FILE* fid = fopen(ss_filename.str().c_str(), "a");
    //cout << "prima for" << endl;

    for(int counter = 0; counter<(int)Lines.size(); counter++)
    {
        // a line
        vecPair vp = Lines[counter];
        Vector2f ex_1 = vp.first;
        Vector2f ex_2 = vp.second;

        // line representation
        Vector4f lineRep = lineRepresentation(vp.first, vp.second);

        // line's polar form
        Vector2f polar = polarRepresentation(vp.first, vp.second);

        // prepare a stringstream
        stringstream ss_lines;

        // write to file
        ss_lines << lineRep(0)<<  "\t" << lineRep(1) << "\t" << lineRep(2) << "\t" << lineRep(3) << "\t" <<
                    polar(0) << "\t" << polar(1) << "\t" <<
                    ex_1(0) << "\t" << ex_1(1) << "\t" <<
                    ex_2(0) << "\t" << ex_2(1) << "\n";
        fputs(ss_lines.str().c_str(), fid);

        if(counter == 0){
            //cout << lineRep<< endl << polar << endl << ex_1 << endl << ex_2 << endl << endl;
        }

    }

}

void printLinesByExtremes(MatrixXf lines, MatrixXf T, string file_name)
{
     // remove and reopen file
    remove(file_name.c_str());
    FILE* fid = fopen(file_name.c_str(), "a");

    // loop over points
    for(int i = 0; i<lines.cols(); i++)
    {
        // create temporary matrices
        MatrixXf temp1 = MatrixXf::Constant(3,1,1);
        MatrixXf temp2 = MatrixXf::Constant(3,1,1);

        // fill temporary matrices
        temp1.block(0,0,2,1) = lines.block(0,i,2,1);
        temp2.block(0,0,2,1) = lines.block(2,i,2,1);

        // products
        Vector3f prod1 = T*temp1;
        Vector3f prod2 = T*temp2;

        // write to file
        stringstream toTheFile;
        toTheFile << prod1(0) << "\t" << prod1(1) << "\n" <<
                     prod2(0) << "\t" << prod2(1) << "\n\n";

        fputs(toTheFile.str().c_str(), fid);

    }

    fclose(fid);
}

void printAssociations(const MatrixXf& Li, const MatrixXf& Lj, const MatrixXf& T, string file_name)
{
    // remove and/or open file
    remove(file_name.c_str());
    FILE* fid = fopen(file_name.c_str(), "a");
    //cout << Lj << endl << endl;

    // loop over cols
    for(int i = 0; i<Li.cols(); i++)
    {
        // get middlepoints
        Vector2f middle_point_i = Li.block(0,i,2,1);
        Vector2f middle_point_j = Lj.block(0,i,2,1);
        //cout << "middlepoints" << endl;
        //cout << middle_point_j << endl << endl;

        // transform middlepoint j
        MatrixXf transf_mp_j = transformVectors(middle_point_j, T);
        //cout << transf_mp_j << endl << endl;

        // prepare output to file
        stringstream toTheFile;
        toTheFile << middle_point_i(0) << "\t" << middle_point_i(1) << "\n" <<
                     transf_mp_j(0,0) << "\t" << transf_mp_j(1,0) << "\n\n";

        // write to file
        fputs(toTheFile.str().c_str(), fid);
    }

    // close file
    fclose(fid);

}

MatrixXf linesByCol(const vector<vecPairsList>& scans_line, int index)
{
    // select a scan's lines
    vecPairsList Lines = scans_line[index];

    // output
    MatrixXf lines_col = MatrixXf::Zero(10,Lines.size());

    for(int counter = 0; counter<(int)Lines.size(); counter++)
    {
        // a line
        vecPair vp = Lines[counter];
        Vector2f ex_1 = vp.first;
        Vector2f ex_2 = vp.second;

        // line representation
        Vector4f lineRep = lineRepresentation(ex_1, ex_2);

        // line's polar form
        Vector2f polar = polarRepresentation(ex_1, ex_2);

        // fill output column
        /*lines_col.block(0,counter,4,1) = lineRep;
        lines_col.block(4,counter,2,1) = polar;
        lines_col.block(6,counter,2,1) = ex_1;
        lines_col.block(8,counter,2,1) = ex_2;*/
        MatrixXf column_ith(10,1);
        column_ith << lineRep(0),lineRep(1), lineRep(2),lineRep(3),
                      polar(0),polar(1),
                      ex_1(0),ex_1(1),
                      ex_2(0),ex_2(1);
        lines_col.block(0,counter,10,1) = column_ith;
    }

    // return output
    return lines_col;
}

MatrixXi findMinIndeces(MatrixXf input_mat)
{
    // prepare output
    MatrixXi output(1,2);

    // find min in input matrix
    float temp_min = input_mat.minCoeff();

    // set a flag
    bool first_row_flag = true;

    // loop over matrix elements
    for(int k = 0; k<input_mat.rows(); k++)
    {
        for (int z = 0; z<input_mat.cols(); z++)
        {
            // if current elemnts holds min values
            if(input_mat(k,z) == temp_min)
            {
                // create a new row
                MatrixXi new_row(1,2);
                new_row << k,z;

                // no row inserted in output matrix yet
                if(first_row_flag)
                {
                    output.row(0) = new_row;
                    first_row_flag = false;
                }

                // at least a row in the matrix
                else
                {
                    MatrixXi temp = MatrixXi::Zero(output.rows()+1,2);
                    temp << output, new_row;
                    output = temp;
                }
            }
        }
    }

    return output;
}

MatrixXi findValueInMatrix(MatrixXi input, float value)
{
    // output
    MatrixXi indeces(1,2);

    // flag
    bool first_row_flag = true;

    for(int i = 0; i<input.rows(); i++)
    {
        for(int j = 0; j<input.cols(); j++)
        {
            if(input(i,j) == value)
            {
                MatrixXi new_row(1,2);
                new_row << i,j;
                if(first_row_flag)
                {
                    indeces.row(0) = new_row;
                    first_row_flag = false;
                }

                else
                {
                    MatrixXi new_matrix(indeces.rows()+1, 2);
                    new_matrix << indeces, new_row;
                    indeces = new_matrix;
                }
            }
        }
    }

    return indeces;
}

void removeRow(MatrixXi& matrix, unsigned int rowToRemove)
{
    unsigned int numRows = matrix.rows()-1;
    unsigned int numCols = matrix.cols();

    if( rowToRemove < numRows )
        matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);

    matrix.conservativeResize(numRows,numCols);
}

MatrixXi findRepetitions(const MatrixXi& assoc)
{
    // temporary matrix
    MatrixXi copy_assoc = MatrixXi::Zero(assoc.rows(),3);
    copy_assoc.block(0,0,assoc.rows(),2) = assoc.block(0,0,assoc.rows(),2);

    for(int i = 0; i<copy_assoc.rows(); i++)
    {
        // occurrences of current element
        if(copy_assoc(i,2) == 0)
        {
            MatrixXi current_indeces = findValueInMatrix(copy_assoc.block(0,1,copy_assoc.rows(),1),copy_assoc(i,1));

            // set flag of occurrences
            if(current_indeces.rows() > 1)
            {
                for(int k = 0; k<current_indeces.rows(); k++)
                    copy_assoc(current_indeces(k,0),2) = 1;
            }
        }
    }

    //cout << copy_assoc << endl << endl;

    // output material
    MatrixXi output(1,2);
    bool first_row_flag = true;

    // delete repetitions in matrix
    for(int j = 0; j<copy_assoc.rows(); j++)
    {
        if(copy_assoc(j,2) == 0)
        {
            MatrixXi new_row(1,2);
            new_row << copy_assoc(j,0), copy_assoc(j,1);

            if(first_row_flag)
            {
                output.row(0) = new_row;
                first_row_flag = false;
            }

            else
            {
                MatrixXi new_matrix(output.rows()+1,2);
                new_matrix << output, new_row;
                output = new_matrix;
            }
        }
    }

    return output;

}


}

#endif // MISC_UTILITIES_H
