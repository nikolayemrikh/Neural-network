//
//  main.cpp
//  neuro_network
//
//  Created by Николай Эмрих on 18.01.17.
//  Copyright © 2017 Николай Эмрих. All rights reserved.
//

#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <random>
#include <cmath>

using namespace std;

static double lowBoundW = -0.3;
static double highBoundW = 0.3;

static double alpha = 0.1;

static size_t n_size = 9;
static size_t m_size = 7;

static vector<vector<int>> number_0 =
{{0, 0, 1, 1, 1, 0, 0},
{0, 1, 0, 0, 0, 1, 0},
{1, 0, 0, 0, 0, 0, 1},
{1, 0, 0, 0, 0, 0, 1},
{1, 0, 0, 0, 0, 0, 1},
{1, 0, 0, 0, 0, 0, 1},
{1, 0, 0, 0, 0, 0, 1},
{0, 1, 0, 0, 0, 1, 0},
{0, 0, 1, 1, 1, 0, 0}};
static vector<double> number_0_correct = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0};

static vector<vector<int>> number_1 =
{{0, 0, 0, 1, 0, 0, 0},
{0, 0, 1, 1, 0, 0, 0},
{0, 1, 0, 1, 0, 0, 0},
{1, 0, 0, 1, 0, 0, 0},
{0, 0, 0, 1, 0, 0, 0},
{0, 0, 0, 1, 0, 0, 0},
{0, 0, 0, 1, 0, 0, 0},
{0, 0, 0, 1, 0, 0, 0},
{0, 0, 0, 1, 0, 0, 0}};
static vector<double> number_1_correct = {0, 1, 0, 0, 0, 0, 0, 0, 0, 0};

static vector<vector<int>> number_2 =
{{0, 0, 1, 1, 1, 0, 0},
{0, 1, 0, 0, 0, 1, 0},
{1, 0, 0, 0, 0, 0, 1},
{1, 0, 0, 0, 0, 1, 0},
{0, 0, 0, 0, 1, 0, 0},
{0, 0, 0, 1, 0, 0, 0},
{0, 0, 1, 0, 0, 0, 0},
{0, 1, 0, 0, 0, 0, 0},
{1, 1, 1, 1, 1, 1, 1}};
static vector<double> number_2_correct = {0, 0, 1, 0, 0, 0, 0, 0, 0, 0};
    
static vector<vector<int>> number_3 =
{{0, 0, 1, 1, 1, 0, 0},
{0, 1, 0, 0, 0, 1, 0},
{1, 0, 0, 0, 0, 0, 1},
{0, 0, 0, 0, 0, 1, 0},
{0, 0, 0, 1, 1, 0, 0},
{0, 0, 0, 0, 0, 1, 0},
{1, 0, 0, 0, 0, 0, 1},
{0, 1, 0, 0, 0, 1, 0},
{0, 0, 1, 1, 1, 0, 0}};
static vector<double> number_3_correct = {0, 0, 0, 1, 0, 0, 0, 0, 0, 0};

static vector<vector<int>> number_4 =
{{0, 0, 0, 0, 1, 0, 1},
{0, 0, 0, 1, 0, 0, 1},
{0, 0, 1, 0, 0, 0, 1},
{0, 1, 0, 0, 0, 1, 1},
{1, 1, 1, 1, 1, 1, 1},
{0, 0, 0, 0, 0, 0, 1},
{0, 0, 0, 0, 0, 0, 1},
{0, 0, 0, 0, 0, 0, 1},
{0, 0, 0, 0, 0, 0, 1}};
static vector<double> number_4_correct = {0, 0, 0, 0, 1, 0, 0, 0, 0, 0};

static vector<vector<int>> number_5 =
{{1, 1, 1, 1, 1, 1, 1},
{1, 0, 0, 0, 0, 0, 0},
{1, 0, 0, 0, 0, 0, 0},
{1, 0, 1, 1, 1, 0, 0},
{1, 1, 0, 0, 0, 1, 0},
{0, 0, 0, 0, 0, 0, 1},
{1, 0, 0, 0, 0, 0, 1},
{0, 1, 0, 0, 0, 1, 0},
{0, 0, 1, 1, 1, 0, 0}};
static vector<double> number_5_correct = {0, 0, 0, 0, 0, 1, 0, 0, 0, 0};

static vector<vector<int>> number_6 =
{{0, 0, 1, 1, 1, 1, 0},
{0, 1, 0, 0, 0, 0, 1},
{1, 0, 0, 0, 0, 0, 0},
{1, 0, 1, 1, 1, 0, 0},
{1, 1, 0, 0, 0, 1, 0},
{1, 0, 0, 0, 0, 0, 1},
{1, 0, 0, 0, 0, 0, 1},
{0, 1, 0, 0, 0, 1, 0},
{0, 0, 1, 1, 1, 0, 0}};
static vector<double> number_6_correct = {0, 0, 0, 0, 0, 0, 1, 0, 0, 0};

static vector<vector<int>> number_7 =
{{1, 1, 1, 1, 1, 1, 1},
{0, 0, 0, 0, 0, 0, 1},
{0, 0, 0, 0, 0, 1, 0},
{0, 0, 0, 0, 1, 0, 0},
{0, 0, 0, 1, 0, 0, 0},
{0, 0, 1, 0, 0, 0, 0},
{0, 1, 0, 0, 0, 0, 0},
{1, 0, 0, 0, 0, 0, 0},
{1, 0, 0, 0, 0, 0, 0}};
static vector<double> number_7_correct = {0, 0, 0, 0, 0, 0, 0, 1, 0, 0};

static vector<vector<int>> number_8 =
{{0, 0, 1, 1, 1, 0, 0},
{0, 1, 0, 0, 0, 1, 0},
{1, 0, 0, 0, 0, 0, 1},
{0, 1, 0, 0, 0, 1, 0},
{0, 0, 0, 1, 0, 0, 0},
{0, 1, 0, 0, 0, 1, 0},
{1, 0, 0, 0, 0, 0, 1},
{0, 1, 0, 0, 0, 1, 0},
{0, 0, 1, 1, 1, 0, 0}};
static vector<double> number_8_correct = {0, 0, 0, 0, 0, 0, 0, 0, 1, 0};

static vector<vector<int>> number_9 =
{{0, 0, 1, 1, 1, 0, 0},
{0, 1, 0, 0, 0, 1, 0},
{1, 0, 0, 0, 0, 0, 1},
{1, 0, 0, 0, 0, 0, 1},
{0, 1, 0, 0, 0, 1, 1},
{0, 0, 1, 1, 1, 0, 1},
{0, 0, 0, 0, 0, 0, 1},
{1, 0, 0, 0, 0, 1, 0},
{0, 1, 1, 1, 1, 0, 0}};
static vector<double> number_9_correct = {0, 0, 0, 0, 0, 0, 0, 0, 0, 1};

static vector<vector<int>> number_3_corrupted =
{{0, 0, 1, 1, 1, 0, 0},
{0, 1, 0, 0, 0, 1, 0},
{1, 0, 0, 0, 0, 0, 1},
{0, 0, 0, 0, 0, 1, 0},
{0, 0, 0, 0, 1, 0, 0},
{0, 0, 1, 1, 0, 0, 0},
{0, 0, 0, 1, 0, 0, 0},
{0, 1, 0, 0, 1, 0, 0},
{0, 0, 1, 1, 0, 0, 0}};

static vector<vector<int>> number_5_corrupted =
{{0, 0, 1, 1, 1, 1, 0},
{0, 1, 0, 0, 0, 0, 1},
{1, 0, 0, 0, 0, 0, 0},
{1, 0, 1, 1, 1, 0, 0},
{1, 1, 0, 0, 0, 1, 0},
{0, 0, 0, 0, 0, 0, 1},
{0, 0, 0, 0, 0, 0, 1},
{0, 0, 0, 0, 0, 1, 0},
{1, 1, 1, 1, 1, 0, 0}};

static vector<vector<int>> number_7_corrupted =
{{1, 1, 1, 1, 1, 1, 1},
{0, 0, 0, 0, 0, 0, 1},
{0, 0, 0, 0, 0, 1, 0},
{0, 0, 0, 0, 1, 0, 0},
{0, 0, 0, 1, 0, 0, 0},
{0, 0, 1, 0, 0, 0, 0},
{0, 1, 0, 0, 0, 0, 0},
{0, 1, 0, 0, 0, 0, 0},
{0, 1, 0, 0, 0, 0, 0}};


class Neuron {
    bool start_;
    bool hidden_;
    bool end_;
public:
    vector<double> W_;
    Neuron() {};
    // Конструктор заполнение весов нейрона
    Neuron(size_t size_W) {
        for (size_t i = 0; i != size_W; ++i) {
            double rd = (double)rand() / RAND_MAX;
            rd = lowBoundW + rd * (highBoundW - lowBoundW);
            W_.push_back(rd);
        };
    }
    
    double get_P(vector<double> & X) {
        double P = 0;
        for (size_t i = 0; i != X.size(); ++i) {
            P += X[i] * W_[i];
        }
        if (P == 0) P = 0;
        return P;
    }
    
    double get_y(double P) {
        double y;
        y = 1 / (1 + exp(-P));
        if (y == 0.0) y = 0.0;
        return y;
    }
};

class NeuroNet {
    vector<Neuron> start_neurons_;
    vector<Neuron> hidden_neurons_;
    vector<Neuron> end_neurons_;
public:
    NeuroNet(size_t start_n, size_t hidden_n, size_t end_n) {
        start_neurons_.resize(start_n);
        hidden_neurons_.resize(hidden_n);
        end_neurons_.resize(end_n);
        for (size_t i = 0; i != hidden_n; ++i) {
            for (size_t j = 0; j != start_n; ++j) {
                hidden_neurons_[i] = Neuron(start_n);
            }
        }
        for (size_t i = 0; i != end_n; ++i) {
            for (size_t j = 0; j != hidden_n; ++j) {
                end_neurons_[i] = Neuron(hidden_n);
            }
        }
    };
    // number - двойной массив представляющий число
    // correct_ans - массив в котором содержатся правильные значения для выходных нейронов
    void teach(vector<vector<int>> & points_array, vector<double> & correct_ans) {
        vector<double> start_y;
        size_t size_s = points_array.size();
        for (size_t i = 0; i != size_s; ++i) {
            size_t size_m = points_array[i].size();
            if (size_m * size_s != start_neurons_.size())
                break;
            for (size_t j = 0; j != size_m; ++j) {
                double val = points_array[i][j];
                start_y.push_back(val);
                //cout << start_y[i * size_m + j];
            }
            //cout << endl;
        }
        vector<double> hidden_y;
        size_t size_hidden = hidden_neurons_.size();
        for (size_t i = 0; i != size_hidden; ++i) {
            double P = hidden_neurons_[i].get_P(start_y);
            double y = hidden_neurons_[i].get_y(P);
            hidden_y.push_back(y);
        }
        
        // Сразу же вычислю ошибку в выходном элементе
        vector<double> end_y;
        vector<double> errors_end;
        size_t size_end = end_neurons_.size();
        for (size_t i = 0; i != size_end; ++i) {
            double P = end_neurons_[i].get_P(hidden_y);
            double y = end_neurons_[i].get_y(P);
            end_y.push_back(y);
            double error = get_error_j(end_y[i], correct_ans[i]);
            errors_end.push_back(error);
        }
        
        // Вычислю ошибку в скрытом слое
        vector<double> errors_hidden;
        for (size_t i = 0; i != size_hidden; ++i) {
            double error = get_error_i(i, hidden_y[i], errors_end, end_neurons_);
            errors_hidden.push_back(error);
        }
        
        // Пересчет весов
        for (size_t i = 0; i != hidden_neurons_.size(); ++i) {
            for (size_t j = 0; j != hidden_neurons_[i].W_.size(); ++j) {
                double val = alpha * errors_hidden[i] * start_y[j];
                if (val == 0.0) val = 0.0;
                hidden_neurons_[i].W_[j] += val;
            }
        }
        
        for (size_t i = 0; i != end_neurons_.size(); ++i) {
            for (size_t j = 0; j != end_neurons_[i].W_.size(); ++j) {
                double val = alpha * errors_end[i] * hidden_y[j];
                if (val == 0.0) val = 0.0;
                end_neurons_[i].W_[j] += val;
            }
        }
        
    }
    
    double get_error_j(double & yj, double & etalon_yj) {
        double val = (etalon_yj - yj) * yj * (1 - yj);
        if (val == 0.0) val = 0.0;
        return val;
    }
    double get_error_i(size_t index, double & yi, vector<double> & errors_j, vector<Neuron> & neurons) {
        size_t size = errors_j.size();
        double sum = 0;
        for (size_t i = 0; i != size; ++i) {
            sum += errors_j[i] * neurons[i].W_[index];
        }
        double val = yi * (1 - yi) * sum;
        if (val == 0.0) val = 0.0;
        return val;
    }
    
    vector<double> get_answer(vector<vector<int>> & points_array) {
        vector<double> start_y;
        size_t size_s = points_array.size();
        for (size_t i = 0; i != size_s; ++i) {
            size_t size_m = points_array[i].size();
            if (size_m * size_s != start_neurons_.size())
                break;
            for (size_t j = 0; j != size_m; ++j) {
                double val = points_array[i][j];
                start_y.push_back(val);
                //cout << start_y[i * size_m + j];
            }
            //cout << endl;
        }
        vector<double> hidden_y;
        size_t size_hidden = hidden_neurons_.size();
        for (size_t i = 0; i != size_hidden; ++i) {
            double P = hidden_neurons_[i].get_P(start_y);
            double y = hidden_neurons_[i].get_y(P);
            hidden_y.push_back(y);
        }
        
        vector<double> end_y;
        size_t size_end = end_neurons_.size();
        for (size_t i = 0; i != size_end; ++i) {
            double P = end_neurons_[i].get_P(hidden_y);
            double y = end_neurons_[i].get_y(P);
            end_y.push_back(y);
        }
        return end_y;
    }
};




// Кол-во нейронов скрытого слоя - 4
// т.е. нужно 25 примеров, содержащих входные данные и правильные ответы
int main() {
    vector<vector<int>> m;
    NeuroNet net(9 * 7, 10, 10);
    for (size_t i = 0; i != 100000; ++i) {
        net.teach(number_0, number_0_correct);
        net.teach(number_1, number_1_correct);
        net.teach(number_2, number_2_correct);
        net.teach(number_3, number_3_correct);
        net.teach(number_4, number_4_correct);
        net.teach(number_5, number_5_correct);
        net.teach(number_6, number_6_correct);
        net.teach(number_7, number_7_correct);
        net.teach(number_8, number_8_correct);
        net.teach(number_9, number_9_correct);
    }
    vector<double> answer3 = net.get_answer(number_3_corrupted);
    vector<double> answer5 = net.get_answer(number_5_corrupted);
    vector<double> answer7 = net.get_answer(number_7_corrupted);
    return 0;
}
