#ifndef MATCHING_UTILITIES_H
#define MATCHING_UTILITIES_H

#include "line_extractor.h"
#include <vector>
#include <map>
#include <Eigen/Core>
#include <Eigen/StdVector>
#include <Eigen/Dense>
#include "stdio.h"
#include <cmath>

#define MAX_LINES_ANGLE_DEG 10
#define MAX_LINES_ANGLE_RAD MAX_LINES_ANGLE_DEG*M_PI/180
#define MAX_LINE_DIST 0.1


using namespace std;
using namespace Eigen;

namespace utilities
{
    typedef map<pair<int, int>, float> linesMap;
    typedef pair<float, float> coordinate;
    typedef pair<Vector2f, Vector2f> vecPair;
    typedef vector<vecPair> vecPairsList;
    struct completeInformation
    {
        coordinate position;
        float angle;
        vector<coordinate> points;
    };

    // compute line representation
    Vector4f lineRepresentation(Vector2f extreme1, Vector2f extreme2)
    {
        // compute tangent vector
        Vector2f tangVec = (extreme2-extreme1)/(extreme2-extreme1).norm();

        // point on the line
        Vector2f pointOnLine = extreme1 + tangVec*(-tangVec.transpose()*extreme1)/tangVec.squaredNorm();

        // compose data
        Vector4f lineRep;
        lineRep << pointOnLine(0), pointOnLine(1), -tangVec(1), tangVec(0);

        // return data
        return lineRep;
    }

    // represent a line in polar form
    Vector2f polarRepresentation(Vector2f extreme1, Vector2f extreme2)
    {
        // local constants
        float x1 = extreme1(0);
        float y1 = extreme1(1);
        float x2 = extreme2(0);
        float y2 = extreme2(1);

        // compute slope
        float angle = atan2(y2-y1, x2-x1);

        // line's coefficients
        float a = y2-y1;
        float b = x1-x2;
        float c = y1*x2 - x1*y2;

        // line's distance from origin
        float rho = (float)fabs(c)/(sqrt(a*a+b*b));

        // output
        Vector2f polar;
        polar << rho, angle;

        return polar;
    }

    // extract lines from a laser scan
    vecPairsList computeLines(Vector2fVector points)
    {
        remove("lines.txt");
        FILE* lines = fopen("lines.txt", "a");
        vecPairsList outputPairs;

        // divide points into clusters
        Point2DClusterer clusterer;
        clusterer.setPoints(points);
        clusterer.compute();

        for (int i = 0; i<clusterer.numClusters(); i++)
        {
            Point2DClusterer::Cluster cluster = clusterer.cluster(i);
            Vector2fVector clusterPoints;
            Line2DExtractor extractor;

            // for current cluster, retrieve its points from input list
            for (int j=cluster.first; j<cluster.second; j++)
            {
                clusterPoints.push_back(points[j]);
            }

            if(clusterPoints.size() > 0)
            {
                // extract lines from point in current cluster
                extractor.setPoints(clusterPoints);
                extractor.compute();

                // save lines extracted from current cluster in outputPairs
                for (Line2DExtractor::IntLineMap::const_iterator it = extractor.lines().begin(); 	 it != extractor.lines().end(); it++)
                {
                     stringstream ss;
                     Vector2f head = clusterPoints[(it->second).p0Index];
                     Vector2f tail = clusterPoints[(it->second).p1Index];

                     outputPairs.push_back(make_pair(head, tail));

                     ss << head(0) << "\t" << head(1) << "\n"
                        << tail(0) << "\t" << tail(1) << "\n\n";

                     fputs(ss.str().c_str(), lines);
                }
            }

        }

        return outputPairs;
    }


}

#endif // MATCHING_UTILITIES_H
