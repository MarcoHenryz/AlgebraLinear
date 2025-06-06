#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
// Definindo nossa própria constante para PI com alta precisão
#define PI 3.14159265358979323846

#define TOLERANCE 1e-9 // Pequena tolerância para comparações de ponto flutuante

// --- ESTRUTURAS E PROTÓTIPOS ---

// Matriz de double para precisão nos cálculos
typedef struct {
    double m[3][3];
} Matriz3x3;

// Vetor de double
typedef struct {
    double v[3];
} Vetor3;

// Protótipos
void imprimir_matriz_int(int **mat1);
void imprimir_matriz_autovalores(double* lambda);
double* polinomio_caracteristico(int **matriz, double* polinomio);
void encontrar_autovalores(double* polinomio, double* lambdas);
void encontrar_autovetores(int **matriz_original, double* lambdas);
void liberar_matriz_int(int **mat);


// --- FUNÇÕES DE IMPRESSÃO ---

void imprimir_matriz_int(int **mat) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%d\t", mat[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void imprimir_matriz_autovalores(double* lambda) {
    printf("Matriz de autovalores (D):\n");
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%.4f\t", (i == j) ? lambda[i] : 0.0);
        }
        printf("\n");
    }
    printf("\n");
}


// --- FUNÇÕES DE CÁLCULO ---

// Calcula os coeficientes do polinômio característico
double* polinomio_caracteristico(int **matriz, double* polinomio) {
    // p(λ) = λ³ - tr(A)λ² + C₂λ - det(A) = 0
    double traco = matriz[0][0] + matriz[1][1] + matriz[2][2];
    double det = (double)matriz[0][0]*(matriz[1][1]*matriz[2][2] - matriz[1][2]*matriz[2][1]) -
                 (double)matriz[0][1]*(matriz[1][0]*matriz[2][2] - matriz[1][2]*matriz[2][0]) +
                 (double)matriz[0][2]*(matriz[1][0]*matriz[2][1] - matriz[1][1]*matriz[2][0]);

    double c2 = (double)(matriz[0][0]*matriz[1][1] - matriz[0][1]*matriz[1][0]) +
                (double)(matriz[0][0]*matriz[2][2] - matriz[0][2]*matriz[2][0]) +
                (double)(matriz[1][1]*matriz[2][2] - matriz[1][2]*matriz[2][1]);
    
    polinomio[0] = 1; // λ³
    polinomio[1] = -traco; // λ²
    polinomio[2] = c2; // λ
    polinomio[3] = -det; // termo constante

    printf("Polinomio caracteristico:\n");
    printf("%.0fx^3 + (%.2f)x^2 + (%.2f)x + (%.2f) = 0\n\n", polinomio[0], polinomio[1], polinomio[2], polinomio[3]);

    return polinomio;
}

// Encontra os autovalores (raízes) usando a solução trigonométrica para cúbicas
void encontrar_autovalores(double* p, double* lambdas) {
    // Normaliza para x³ + ax² + bx + c = 0
    double a = p[1] / p[0];
    double b = p[2] / p[0];
    double c = p[3] / p[0];

    // Converte para a forma y³ + py + q = 0
    double q_form = (3*b - a*a) / 9.0;
    double r_form = (9*a*b - 27*c - 2*a*a*a) / 54.0;
    double discriminante = q_form*q_form*q_form + r_form*r_form;

    if (discriminante < 0) { // Três raízes reais (caso comum para matrizes simétricas)
        double theta = acos(r_form / sqrt(-(q_form*q_form*q_form)));
        double sqrt_q = 2 * sqrt(-q_form);
        // Usa a nossa constante PI definida manualmente no lugar de M_PI
        lambdas[0] = sqrt_q * cos(theta / 3.0) - a / 3.0;
        lambdas[1] = sqrt_q * cos((theta + 2.0*PI) / 3.0) - a / 3.0;
        lambdas[2] = sqrt_q * cos((theta + 4.0*PI) / 3.0) - a / 3.0;
    } else { // Uma raiz real e duas complexas (ou raízes reais repetidas)
        double S = cbrt(r_form + sqrt(discriminante));
        double T = cbrt(r_form - sqrt(discriminante));
        lambdas[0] = (S + T) - a / 3.0;
        lambdas[1] = -(S + T)/2.0 - a / 3.0; 
        lambdas[2] = lambdas[1]; // Simplificação: mostra apenas a parte real
    }
    
    printf("Autovalores encontrados:\n");
    printf("lambda[1] = %.4f\n", lambdas[0]);
    printf("lambda[2] = %.4f\n", lambdas[1]);
    printf("lambda[3] = %.4f\n\n", lambdas[2]);
}

// Encontra os autovetores associados a cada autovalor
void encontrar_autovetores(int **matriz_original, double* lambdas) {
    printf("Autovetores associados a cada autovalor:\n");

    for (int i = 0; i < 3; i++) {
        double lambda = lambdas[i];
        Matriz3x3 B; // Matriz (A - λI)
        
        // Monta a matriz B
        for (int row = 0; row < 3; row++) {
            for (int col = 0; col < 3; col++) {
                B.m[row][col] = matriz_original[row][col];
            }
            B.m[row][row] -= lambda;
        }

        // Calcula o produto vetorial de duas linhas da matriz B
        // para encontrar um vetor no espaço nulo (o autovetor).
        Vetor3 v1 = {B.m[0][0], B.m[0][1], B.m[0][2]};
        Vetor3 v2 = {B.m[1][0], B.m[1][1], B.m[1][2]};
        Vetor3 autovetor;

        autovetor.v[0] = v1.v[1]*v2.v[2] - v1.v[2]*v2.v[1];
        autovetor.v[1] = v1.v[2]*v2.v[0] - v1.v[0]*v2.v[2];
        autovetor.v[2] = v1.v[0]*v2.v[1] - v1.v[1]*v2.v[0];
        
        // Normaliza o autovetor (torna seu comprimento igual a 1)
        double mag = sqrt(autovetor.v[0]*autovetor.v[0] + autovetor.v[1]*autovetor.v[1] + autovetor.v[2]*autovetor.v[2]);

        if (mag > TOLERANCE) {
            autovetor.v[0] /= mag;
            autovetor.v[1] /= mag;
            autovetor.v[2] /= mag;
        }

        printf("Autovetor para lambda = %.4f: [%.4f, %.4f, %.4f]\n", lambda, autovetor.v[0], autovetor.v[1], autovetor.v[2]);
    }
    printf("\n");
}

void liberar_matriz_int(int **mat) {
    if (mat == NULL) return;
    for (int i = 0; i < 3; i++) {
        free(mat[i]);
    }
    free(mat);
}


// --- FUNÇÃO PRINCIPAL ---
int main(int argc, char** argv) {
    if (argc < 2) {
        fprintf(stderr, "Erro: Forneca o nome do arquivo de entrada como argumento.\n");
        fprintf(stderr, "Uso: %s <arquivo_de_entrada.txt>\n", argv[0]);
        return 1;
    }

    FILE* input = fopen(argv[1], "r");
    if (input == NULL) {
        fprintf(stderr, "Erro ao abrir o arquivo de entrada '%s'\n", argv[1]);
        return 1;
    }

    if (freopen("danilokotakaemarcohenry.txt", "w", stdout) == NULL) {
        fprintf(stderr, "Erro ao criar o arquivo de saida.\n");
        fclose(input);
        return 1;
    }

    printf("Nomes: Danilo Kotaka e Marco Henry\n\n");

    int** matriz = calloc(3, sizeof(int*));
    for (int i = 0; i < 3; i++) {
        matriz[i] = calloc(3, sizeof(int));
        for (int j = 0; j < 3; j++) {
            fscanf(input, "%d", &matriz[i][j]);
        }
    }
    fclose(input);

    printf("Matriz de entrada (A):\n");
    imprimir_matriz_int(matriz);

    double polinomio[4];
    double lambdas[3];
    
    polinomio_caracteristico(matriz, polinomio);
    encontrar_autovalores(polinomio, lambdas);
    imprimir_matriz_autovalores(lambdas);
    encontrar_autovetores(matriz, lambdas);

    liberar_matriz_int(matriz);

    return 0;
}