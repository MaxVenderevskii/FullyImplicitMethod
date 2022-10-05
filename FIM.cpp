
#include <iostream>
#include <vector>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

double vol(double a, double b) {
	return a * b;
}

double sigma(double x, double y) {
	return x / y / y;
}

double harm_mean_k(double x, double y) {
	return x*y/(x+y);
}


struct system_status
{
    double tau;
    double hx;
    double hy;
    int M; // Ðàçìåðíîñòü ïî x
    int N; // Ðàçìåðíîñòü ïî y
    double boundary_condition_left; // Çíà÷åíèå äëÿ ãðàíè÷íîãî óñëîâèÿ I ðîäà ñëåâà
    double boundary_condition_right; // Çíà÷åíèå äëÿ ãðàíè÷íîãî óñëîâèÿ I ðîäà ñïðàâà
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

    std::vector<std::vector<double>> matrix(M * N, std::vector<double>(M * N, 0)); // Ââîäèì ìàòðèöó êîýôôèöèåíòîâ äëÿ âñåõ òî÷åê ñåòêè

    std::vector<double> right_parts(M * N, 0); // Ââîäèì âåêòîð, ñîäåðæàùèé ãðàíè÷íûå óñëîâèÿ, äàëåå, ïðè íåîáõîäèìîñòè çàïîëíÿåì çíàåíèÿìè u[i, j] íà ïðåäûäóùåì øàãå

    auto ind = [&](int i, int j) {
        return j * M + i;
    };

    // Çàïîëíÿåì ìàòðèöó MNxMN
    for (int I = 0; I < M * N; I++) {

        int i = I % M;
        int j = I / M;

        matrix[I][I] = 1 / tau;


        if (i > 0) { // поток слева направо
            matrix[I][ind(i - 1, j)] = ((-2) * k_mean(k_matrix[j][i],k_matrix[j][i - 1]) / (hx * hx); // [I - 1] 
            matrix[I][I] += (k_mean(k_matrix[j][i],k_matrix[j][i - 1]) / (hx * hx);
        }
        if (i < M - 1) { // поток справа налево
            matrix[I][ind(i + 1, j)] = ((-2) * k_mean(k_matrix[j][i], k_matrix[j][i + 1]) / (hx * hx); // [I + 1]
            matrix[I][I] += (k_mean(k_matrix[j][i],k_matrix[j][i + 1]) / (hx * hx);
        }
        if (j > 0) { // поток снизу вверх
            matrix[I][ind(i, j - 1)] = ((-2) * k_mean(k_matrix[j - 1][i], k_matrix[j][i]) / (hy * hy);
            matrix[I][I] += (k_mean(k_matrix[j][i],k_matrix[j - 1][i]) / (hy * hy);
        }
        if (j < N - 1) { // поток сверху вниз
            matrix[I][ind(i, j + 1)] = ((-2) * k_mean(k_matrix[j + 1][i], k_matrix[j][i]) / (hy * hy);
            matrix[I][I] += (k_mean(k_matrix[j][i], k_matrix[j + 1][i]) / (hy * hy);
        }

        if (i == 0) { // граничные условия слева
            right_parts[I] += bond_cond_l / (hx * hx);
            matrix[I][I] += 2*k_matrix[i][j]/(hx * hx);  
        }
        if (i == M - 1) { // граничные условия справа
            right_parts[I] += bond_cond_r / (hx * hx);
            matrix[I][I] 
        }

    }

    // âûâîä ìàññèâà
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

    // çàïîëíÿåì ìàòðèöó êîýôôèöèåíòîâ k
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

// Çàïóñê ïðîãðàììû: CTRL+F5 èëè ìåíþ "Îòëàäêà" > "Çàïóñê áåç îòëàäêè"
// Îòëàäêà ïðîãðàììû: F5 èëè ìåíþ "Îòëàäêà" > "Çàïóñòèòü îòëàäêó"

// Ñîâåòû ïî íà÷àëó ðàáîòû 
//   1. Â îêíå îáîçðåâàòåëÿ ðåøåíèé ìîæíî äîáàâëÿòü ôàéëû è óïðàâëÿòü èìè.
//   2. Â îêíå Team Explorer ìîæíî ïîäêëþ÷èòüñÿ ê ñèñòåìå óïðàâëåíèÿ âåðñèÿìè.
//   3. Â îêíå "Âûõîäíûå äàííûå" ìîæíî ïðîñìàòðèâàòü âûõîäíûå äàííûå ñáîðêè è äðóãèå ñîîáùåíèÿ.
//   4. Â îêíå "Ñïèñîê îøèáîê" ìîæíî ïðîñìàòðèâàòü îøèáêè.
//   5. Ïîñëåäîâàòåëüíî âûáåðèòå ïóíêòû ìåíþ "Ïðîåêò" > "Äîáàâèòü íîâûé ýëåìåíò", ÷òîáû ñîçäàòü ôàéëû êîäà, èëè "Ïðîåêò" > "Äîáàâèòü ñóùåñòâóþùèé ýëåìåíò", ÷òîáû äîáàâèòü â ïðîåêò ñóùåñòâóþùèå ôàéëû êîäà.
//   6. ×òîáû ñíîâà îòêðûòü ýòîò ïðîåêò ïîçæå, âûáåðèòå ïóíêòû ìåíþ "Ôàéë" > "Îòêðûòü" > "Ïðîåêò" è âûáåðèòå SLN-ôàéë.
