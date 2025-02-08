#include <bits/stdc++.h>
#include <glad/glad.h>
#include <matplot/matplot.h>
#include <matplot/backend/opengl.h>
#include <GLFW/glfw3.h>

using namespace std;
using namespace matplot;

int GROUP = 4;
int X = 4;
int Y = 4;
int Z = 4;
int N = GROUP * X * Y * Z;
double A = 2;

double initial[4][3] = {
    {0, 0, 0}, {A / 2, A / 2, 0},
    {0, A / 2, A / 2}, {A / 2, 0, A / 2}
};

struct molecule {
    int property;
    double x, y, z;
    double vx, vy, vz;
};

int main() {
    // GLFWの初期化
    if (!glfwInit()) {
        cerr << "Failed to initialize GLFW" << endl;
        return -1;
    }
/*
    // OpenGLコンテキストのバージョンを設定（4.4）
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 4);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#failed to create window
g++ test.cpp -o test -std=c++17 \
    -I/home/amon/vcpkg/installed/x64-linux/include \
    -L/home/amon/vcpkg/installed/x64-linux/lib \
    -lmatplot -lnodesoup -lmatplot_opengl \
    -lglfw3 -lglad -lGL -ldl -lstdc++fs -lpthread && ./test*/

    glfwWindowHint(GLFW_RESIZABLE, GLFW_FALSE);

    // ウィンドウ作成
    GLFWwindow* window = glfwCreateWindow(800, 600, "Test Window", NULL, NULL);
    if (!window) {
        const char* description;
        int code = glfwGetError(&description);
        cerr << "Failed to create GLFW window. Error Code: " << code
             << ", Description: " << (description ? description : "Unknown") << endl;
        glfwTerminate();
        return -1;
    }


    // OpenGLのコンテキストを設定
    glfwMakeContextCurrent(window);
    if (!gladLoadGLLoader(reinterpret_cast<GLADloadproc>(glfwGetProcAddress))) {
        cerr << "Failed to initialize GLAD" << endl;
        return -1;
    }

    // Matplot++のOpenGLバックエンドを使用
    auto f = figure(true);
    f->backend(std::make_shared<backend::opengl>());

    // 分子の配置
    vector<molecule> fcc(N);
    for (int i = 0; i < N / GROUP; i++) {
        for (int j = 0; j < GROUP; j++) {
            fcc[i * GROUP + j].property = j;
            fcc[i * GROUP + j].x = initial[j][0] + (i % X) * A;
            fcc[i * GROUP + j].y = initial[j][1] + ((i / X) % Y) * A;
            fcc[i * GROUP + j].z = initial[j][2] + ((i / (X * Y)) % Z) * A;
        }
    }

    // 速度の初期化（全体の運動量をゼロにする）
    mt19937_64 mt64(0);
    uniform_real_distribution<double> uni(-1, 1);
    double total_vx = 0, total_vy = 0, total_vz = 0;

    for (int i = 0; i < N; i++) {
        fcc[i].vx = uni(mt64);
        fcc[i].vy = uni(mt64);
        fcc[i].vz = uni(mt64);
        total_vx += fcc[i].vx;
        total_vy += fcc[i].vy;
        total_vz += fcc[i].vz;
    }

    for (int i = 0; i < N; i++) {
        fcc[i].vx -= total_vx / N;
        fcc[i].vy -= total_vy / N;
        fcc[i].vz -= total_vz / N;
    }

    // OpenGLの描画ループ
    while (!glfwWindowShouldClose(window)) {
        // 背景をクリア
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // 分子の座標を取得
        vector<double> x(N), y(N), z(N);
        transform(fcc.begin(), fcc.end(), x.begin(), [](const molecule& p) { return p.x; });
        transform(fcc.begin(), fcc.end(), y.begin(), [](const molecule& p) { return p.y; });
        transform(fcc.begin(), fcc.end(), z.begin(), [](const molecule& p) { return p.z; });

        // Matplot++を使用して描画
        scatter3(x, y, z);

        // OpenGLのバッファを交換
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    // GLFWの終了処理
    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}
