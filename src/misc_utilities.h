#ifndef MISC_UTILITIES_H
#define MISC_UTILITIES_H

#include "matching_utilities.h"
#include "fstream"

#define MIN_LINE_LENGTH 0.1*0.1
#define MIN_FRAME_DIST 0.1*0.1

namespace utilities
{

void parseRobotInfo(char* file_name, vector<completeInformation>& ci)
{
    ifstream file(file_name);
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
        cout << "Remaining lines: " << new_vpl.size() << "/" << vpl.size() << endl;
    }

    return selectedLines;
}

void printLines(vecPairsList Lines, int index)
{
    cout << "entrata" << endl;
    stringstream ss_filename;
    ss_filename << "convertedLines_" << index << ".txt";

    remove(ss_filename.str().c_str());
    FILE* fid = fopen(ss_filename.str().c_str(), "a");
    cout << "prima for" << endl;

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

    }

}

}

#endif // MISC_UTILITIES_H
