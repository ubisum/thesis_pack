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

Vector4f lsIteration(MatrixXf X, MatrixXf Z, float epsilon, float alpha)
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
        MatrixXf r_i(2,2);
        r_i << normal(0), -normal(1), normal(1), normal(0);

        // covariance matrices
        MatrixXf epsilon_matrix = MatrixXf::Identity(2,2);
        epsilon_matrix(1,1) = epsilon;

        MatrixXf sigma_pi = r_i*epsilon_matrix*r_i.transpose();

        MatrixXf sigma_ij = MatrixXf::Zero(4,4);
        sigma_ij.block(0,0,2,2) = sigma_pi;
        sigma_ij.block(2,2,2,2) = alpha*MatrixXf::Identity(2,2);

        // local matrices
        MatrixXf h_ij = jacobian.transpose()*sigma_ij*jacobian;
        Vector3f b_ij = jacobian.transpose()*sigma_ij*error;

        // update global matrices
        H+= h_ij;
        b+=b_ij;
        chi+= error.transpose()*error;
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

}

#endif // LEAST_SQUARES_UTILITIES_H
