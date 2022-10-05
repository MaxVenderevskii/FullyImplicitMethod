
#include <iostream>
#include <vector>
#include <Eigen/Dense>


struct system_status
{
    double tau;
    double hx;
    double hy;
    int M; // ����������� �� x
    int N; // ����������� �� y
    double boundary_condition_left; // �������� ��� ���������� ������� I ���� �����
    double boundary_condition_right; // �������� ��� ���������� ������� I ���� ������
};

std::pair<std::vector<std::vector<double>>, std::vector<double>> matrix_equation(const std::vector<std::vector<double>>& k_matrix, system_status constant, bool bond_cond)
{
    double tau = constant.tau;
    double hx = constant.hx;
    double hy = constant.hy;
    int M = constant.M;
    int N = constant.N;
    double bond_cond_l = constant.boundary_condition_left;
    double bond_cond_r = constant.boundary_condition_right;

    std::vector<std::vector<double>> matrix(M * N, std::vector<double>(M * N, 0)); // ������ ������� ������������� ��� ���� ����� �����

    std::vector<double> right_parts(M * N, 0); // ������ ������, ���������� ��������� �������, �����, ��� ������������� ��������� ��������� u[i, j] �� ���������� ����

    int i = 0, j = 0;

    auto ind = [&](int i, int j) {
        return j * M + i;
    };

    // ��������� ������� MNxMN
    for (int I = 0; I < M * N; I++) {

        i = I % M;
        j = I / M;

        matrix[I][I] = 1 / tau;


        if (i > 0) {
            matrix[I][ind(i - 1, j)] = ((-2) / (1 / k_matrix[j][i] + 1 / k_matrix[j][i - 1])) / (hx * hx); // [I - 1] 
            matrix[I][I] += (1 / (1 / k_matrix[j][i] + 1 / k_matrix[j][i - 1])) / (hx * hx);
        }
        if (i < M - 1) {
            matrix[I][ind(i + 1, j)] = ((-2) / (1 / k_matrix[j][i] + 1 / k_matrix[j][i + 1])) / (hx * hx); // [I + 1]
            matrix[I][I] += (1 / (1 / k_matrix[j][i] + 1 / k_matrix[j][i + 1])) / (hx * hx);
        }
        if (j > 0) {
            matrix[I][ind(i, j - 1)] = ((-2) / (1 / k_matrix[j - 1][i] + 1 / k_matrix[j][i])) / (hy * hy);
            matrix[I][I] += (1 / (1 / k_matrix[j][i] + 1 / k_matrix[j - 1][i])) / (hy * hy);
        }
        if (j < N - 1) {
            matrix[I][ind(i, j + 1)] = ((-2) / (1 / k_matrix[j + 1][i] + 1 / k_matrix[j][i])) / (hy * hy);
            matrix[I][I] += (1 / (1 / k_matrix[j][i] + 1 / k_matrix[j + 1][i])) / (hy * hy);
        }

        if (i == 0) {
            right_parts[I] += bond_cond_l / (hx * hx);
        }
        if (i == M - 1) {
            right_parts[I] += bond_cond_r / (hx * hx);
        }

    }

    // ����� �������
    for (int p = 0; p < M * N; p++) {
        for (int q = 0; q < M * N; q++) {
            std::cout << matrix[p][q] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;

    for (int I = 0; I < M * N; I++) {
        std::cout << right_parts[I] << std::endl;
    }


    return std::make_pair(matrix, right_parts);


}




int main()
{
    system_status constant = { 0.1, 0.2, 0.3, 3, 3, 0., 0. }; // {tau, hx, hy, M, N, bond_cond_l, bond_cond_r}

    bool boundary_conditions_I = true;

    int M = constant.M;
    int N = constant.N;

    // ��������� ������� ������������� k
    std::vector<std::vector<double>> k_matrix(N, std::vector<double>(M, 1));
    for (int p = 0; p < N; p++) {
        for (int q = 0; q < M; q++) {
            k_matrix[p][q] = (q < (M / 2 + 1) ? 2 : 3);
        }
    }

    auto equation = matrix_equation(k_matrix, constant, boundary_conditions_I);

}





void func_matrix() {
    Eigen::Matrix3f A;
    Eigen::Vector3f b;
    A << 1, 2, 3, 4, 5, 6, 7, 8, 10;
    b << 3, 3, 4;
    std::cout << "Here is the matrix A:\n" << A << std::endl;
    std::cout << "Here is the vector b:\n" << b << std::endl;
    Eigen::Vector3f x = A.colPivHouseholderQr().solve(b);
    std::cout << "The solution is:\n" << x << std::endl;
}

// ������ ���������: CTRL+F5 ��� ���� "�������" > "������ ��� �������"
// ������� ���������: F5 ��� ���� "�������" > "��������� �������"

// ������ �� ������ ������ 
//   1. � ���� ������������ ������� ����� ��������� ����� � ��������� ���.
//   2. � ���� Team Explorer ����� ������������ � ������� ���������� ��������.
//   3. � ���� "�������� ������" ����� ������������� �������� ������ ������ � ������ ���������.
//   4. � ���� "������ ������" ����� ������������� ������.
//   5. ��������������� �������� ������ ���� "������" > "�������� ����� �������", ����� ������� ����� ����, ��� "������" > "�������� ������������ �������", ����� �������� � ������ ������������ ����� ����.
//   6. ����� ����� ������� ���� ������ �����, �������� ������ ���� "����" > "�������" > "������" � �������� SLN-����.
