#ifndef LEAST_SQUARES_UTILITIES_H
#define LEAST_SQUARES_UTILITIES_H



namespace least_squares
{

/* X is a 3x3 transformation matrix, while P is a 4xN matrix, containing
lines by columns, in the form [point, normal]' */
MatrixXf transformPoints(MatrixXf X, MatrixXf P)
{
    // create a new matrix to host points
    MatrixXf tempMatrix = MatrixXf::Zero(4, P.cols());

    // extract rotational and traslational parts from X
    MatrixXf rotmat = X.block<2,2>(0,0);
    Vector2f traMat = X.block<2,1>(0,2);

    // separate points from normals
    MatrixXf points = P.block(0,0,2,P.cols());
    MatrixXf normals = P.block(2,0,2,P.cols());

    // fill temp matrix
    tempMatrix.block(0,0,2,P.cols()) = rotmat*points + traMat.replicate(1,P.cols());
    tempMatrix.block(2,0,2,P.cols()) = rotmat*normals;

    // return transformed points
    return tempMatrix;
}

/* X is a 3x3 transformation matrix, while Z is a 8xN matrix, containing
couples of lines by columns, each line in the form [point, normal]' */
Vector4f computeError(int index, MatrixXf X, MatrixXf Z)
{
    // get first line
    Vector4f pi = Z.block(0, index, 4, 1);

    // get second line
    Vector4f pj = Z.block(4, index, 4, 1);

    // return error
    return transformPoints(X, pj) - pi;
}

/* X is a 3x3 transformation matrix, while Z is a 8xN matrix, containing
couples of lines by columns, each line in the form [point, normal]' */
MatrixXf computeJacobian (int index, MatrixXf X, MatrixXf Z)
{
    // prepare jacobian
    MatrixXf jacobian = MatrixXf::Zero(4,3);
    jacobian.block(0,0,2,2) = MatrixXf::Identity(2,2);

    // extract i-th line
    Vector4f line_ith = Z.block(4,index,4,1);

    // transform line
    Vector4f tr_points = transformPoints(X, line_ith);

    // separate parts in transformed line
    Vector2f p_i = tr_points.block(0,0,2,1);
    Vector2f n_i = tr_points.block(2,0,2,1);

    // fill last column of Jacobian
    Vector4f last_col(4,1);
    last_col << -p_i(1), p_i(0), -n_i(1), n_i(0);
    jacobian.block(0,2,4,1) = last_col;

    // return Jacobian
    return jacobian;


}

Vector4f lsIteration(MatrixXf X, MatrixXf Z, float epsilon, float alpha, float gamma=1)
{
    // initialize matrices and variables
    MatrixXf H = MatrixXf::Zero(3,3);
    Vector3f b = Vector3f::Zero(3,1);
    float chi = 0;

    // loop
    for(int i = 0; i<Z.cols(); i++)
    {
        // compute error
        Vector4f error = computeError(i,X,Z);

        // compute jacobian
        MatrixXf jacobian = computeJacobian(i,X,Z);

        // compute matrix for line rotation
        Vector2f normal = Z.block(2,i,2,1);
        //normal /= sqrt(normal.transpose()*normal);
        float normal_factor = normal.transpose()*normal;
        normal /= normal_factor;
        MatrixXf r_i(2,2);
        r_i << normal(0), -normal(1), normal(1), normal(0);

        // covariance matrices
        MatrixXf epsilon_matrix = MatrixXf::Identity(2,2);
        epsilon_matrix(1,1) = epsilon;
        epsilon_matrix(0,0) = 10;

        MatrixXf sigma_pi = r_i.transpose()*epsilon_matrix*r_i;

        MatrixXf sigma_ij = MatrixXf::Zero(4,4);
        sigma_ij.block(0,0,2,2) = sigma_pi;
        sigma_ij.block(2,2,2,2) = alpha*MatrixXf::Identity(2,2);

        // scale covariance matrix
        float point_error = error.transpose()*sigma_ij*error;
        float scale = 1;
        if(point_error > gamma)
            scale = sqrt(gamma/point_error);
        sigma_ij *= scale;

        // local matrices
        MatrixXf h_ij = jacobian.transpose()*sigma_ij*jacobian;
        Vector3f b_ij = jacobian.transpose()*sigma_ij*error;

        // update global matrices
        H += h_ij;
        b += b_ij;
        chi += error.transpose()*sigma_ij*error;
    }

    // solve linear system
    //cout << H << endl << endl;
    //cout << b << endl << endl;
    Vector3f dx = (-H).colPivHouseholderQr().solve(b);
    cout << dx << endl << endl;
    MatrixXf DX(3,3);
    DX << cos(dx(2)), -sin(dx(2)), dx(0),
            sin(dx(2)), cos(dx(2)), dx(1),
            0, 0, 1;

    // compute new transformation matrix
    MatrixXf new_transf_mat = DX*X;

    // prepare output
    Vector4f output(4,1);
    output << new_transf_mat(0,2), new_transf_mat(1,2),
            atan2(new_transf_mat(1,0), new_transf_mat(0,0)), chi;

    // return output
    return output;

}

MatrixXf computeDistanceNM(MatrixXf ne_i, MatrixXf ne_j, float alpha_factor, float point_factor, MatrixXf T)
{
    // prepare output
    MatrixXf dist = MatrixXf::Zero(ne_i.cols(), ne_j.cols());

    // transform ne_j
    MatrixXf transf = T.block(0,0,2,2)*ne_j.block(0,0,2,ne_j.cols());
    MatrixXf transf_ex1 = T.block(0,0,2,2)*ne_j.block(2,0,2,ne_j.cols()) +
            (T.block(0,2,2,1)).replicate(1,ne_j.cols());
    MatrixXf transf_ex2 = T.block(0,0,2,2)*ne_j.block(4,0,2,ne_j.cols()) +
            (T.block(0,2,2,1)).replicate(1,ne_j.cols());

    // construct new ne_j
    MatrixXf ne_j_transf = MatrixXf::Zero(6,ne_j.cols());
    ne_j_transf.block(0,0,2,ne_j.cols()) = transf;
    ne_j_transf.block(2,0,2,ne_j.cols()) = transf_ex1;
    ne_j_transf.block(4,0,2,ne_j.cols()) = transf_ex2;

    // loop over lines in ne_i
    for(int i = 0; i<ne_i.cols(); i++)
    {
        // get a column of a line
        MatrixXf col_i = ne_i.block(0,i,6,1); // 6x1

        // compute middle point
        Vector2f mid_i;
        mid_i << col_i(2,0)+col_i(4,0), col_i(3,0)+col_i(5,0);
        mid_i *= 0.5;

        // normal part
        Vector2f n_i = ne_i.block(0,i,2,1);

        // loop over lines in ne_j
        for(int j = 0; j<ne_j.cols(); j++)
        {
            // get a column of line
            MatrixXf col_j = ne_j.block(0,j,6,1); // 6x1

            // compute middle point
            Vector2f mid_j;
            mid_j << col_j(2,0)+col_j(4,0), col_j(3,0)+col_j(5,0);
            mid_j *= 0.5;

            // normal part
            Vector2f n_j = ne_j.block(0,j,2,1);

            // distances
            Vector2f dist_n = n_i-n_j;
            Vector2f dist_m = mid_i-mid_j;
            float dist_n_value = (dist_n.transpose()*dist_n);
            float dist_m_value = (dist_m.transpose()*dist_m);

            // set value
            dist(i,j) = alpha_factor*dist_n_value + point_factor*dist_m_value;
        }


    }

    // return matrix of distances
    return dist;
}

}

#endif // LEAST_SQUARES_UTILITIES_H
