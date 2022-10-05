
#include <iostream>
#include <vector>
#include <Eigen/Dense>


struct system_status
{
    double tau;
    double hx;
    double hy;
    int M; // Размерность по x
    int N; // Размерность по y
    double boundary_condition_left; // Значение для граничного условия I рода слева
    double boundary_condition_right; // Значение для граничного условия I рода справа
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

    std::vector<std::vector<double>> matrix(M * N, std::vector<double>(M * N, 0)); // Вводим матрицу коэффициентов для всех точек сетки

    std::vector<double> right_parts(M * N, 0); // Вводим вектор, содержащий граничные условия, далее, при необходимости заполняем знаениями u[i, j] на предыдущем шаге

    int i = 0, j = 0;

    auto ind = [&](int i, int j) {
        return j * M + i;
    };

    // Заполняем матрицу MNxMN
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

    // вывод массива
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

    // заполняем матрицу коэффициентов k
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

// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
