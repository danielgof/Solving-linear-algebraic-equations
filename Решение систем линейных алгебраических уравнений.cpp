#include <iostream>
#include <cmath>
#define EPS 0.000001

void SwapLine(double* mat, int size, int l1, int l2)
{
    for (int j = 0; j < size; j++)
    {
        double temp = mat[l1 * size + j];
        mat[l1 * size + j] = mat[l2 * size + j];
        mat[l2 * size + j] = temp;
    }
}


void SubMatr(double* mat, int size, int l1, double k1, int l2, double k2, int s)
{

    for (int j = s; j < size; j++)
    {
        mat[l1 * size + j] = mat[l1 * size + j] * k2 - mat[l2 * size + j] * k1;
    }
}

double Gauss(double* mat, int row, int col)
{
    int size = col + 1;

    double koef = 1;
    if (fabs(mat[0 * size + 0]) < EPS)
    {
        for (int i = 1; i < row; i++)
        {
            if (fabs(mat[i * size + 0]) > EPS)
            {
                SwapLine(mat, size, 0, i);
                koef *= -1.0;
                break;
            }
        }
    }

    for (int i0 = 1; i0 < col; i0++)
    {
        for (int ll = i0 + 1; ; ll++)
        {
            for (int j = 0; j < i0; j++)
            {
                if (fabs(mat[i0 * size + j]) < EPS)
                    continue;

                koef *= mat[j * size + j];
                SubMatr(mat, size, i0, mat[i0 * size + j], j, mat[j * size + j], j);
            }
            if (fabs(mat[i0 * size + i0]) > EPS)
            {

                break;
            }
            if (ll == row)
            {
                return 0;
            }

            SwapLine(mat, size, i0, ll);

            koef *= -1;
        }
    }



    for (int i = col - 1; i > 0; i--)
    {
        for (int j = i - 1; j >= 0; j--)
        {
            // mat[j] - mat[i]
            SubMatr(mat, size, j, mat[j * size + i], i, mat[i * size + i], j);
            //SubMatr(mat, size, i0, mat[i0 * size + j], j, mat[j * size + j], j);
        }
    }

    return 1;

}




int Search(double* m, int col, int row, int* uI, int* uJ, int s)
{
    for (int i = s; i < row; i++)
    {
        for (int j = s; j < col; j++)
        {
            if (m[i * col + j] != 0)
            {
                *uI = i;
                *uJ = j;
                return 0;
            }
        }
    }

    return 1;
}

int Rang(double* m, int row, int col)
{
    if (col == 1)
    {
        return 1;
    }

    double tmp;
    int i;
    int x, y;
    // выбираем минимум
    int endi = (col > row) ? row : col;
    // цикл по всей главной диагонали
    for (i = 0; i < endi; i++)
    {
        // если элемент на диагонали равен 0, то ищем не нулевой элемент в матрице
        if (m[i * col + i] == 0)
        {
            // если все элементы нулевые, прерываем цикл
            if (Search(m, col, row, &y, &x, i))
                break;


            if (i != y)
            {
                for (int j = 0; j < col; j++)
                {
                    double temp = m[i * col + j];
                    m[i * col + j] = m[y * col + j];
                    m[y * col + j] = temp;
                }
                //m.swaprows(i, y);
            }
            // меняем i-ый столбец с x-ым
            if (i != x)
            {
                for (int j = 0; j < row; j++)
                {
                    double temp = m[j * col + i];
                    m[j * col + i] = m[j * col + x];
                    m[j * col + x] = temp;
                }
                //m.swapcolumns(i, x);
            }
            // таким образом, в m[i][i], теперь ненулевой элемент.
        }

        // выносим элемент m[i][i]
        tmp = m[i * col + i];
        for (x = i; x < col; x++)
        {
            m[i * col + x] = m[i * col + x] / tmp;
        }
        // таким образом m[i][i] теперь равен 1
        // зануляем все элементы стоящие под (i, i)-ым и справа от него,
        // при помощи вычитания с опр. коеффициентом
        for (y = i + 1; y < row; y++)
        {
            tmp = m[y * col + i];
            for (x = i; x < col; x++)
                m[y * col + x] -= (m[i * col + x] * tmp);
        }
        for (x = i + 1; x < col; x++)
        {
            tmp = m[i * col + x];
            for (y = i; y < row; y++)
                m[y * col + x] -= (m[y * col + i] * tmp);
        }
    }

    // считаем сколько единичек на главной диагонали
    unsigned cnt = 0;
    for (i = 0; i < endi; i++)
        if (m[i * col + i] == 0)
            break;
        else cnt++;

    if (!cnt)
        cnt++;

    return cnt;
}

int main()
{
    int row, col;
    std::cin >> row >> col;

    double* full = new double[row * (col + 1)];
    double* full_cpy = new double[row * (col + 1)];
    double* mat = new double[row * col];


    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            std::cin >> mat[i * col + j];
            full[i * (col + 1) + j] = mat[i * col + j];
            full_cpy[i * (col + 1) + j] = full[i * (col + 1) + j];
        }

        std::cin >> full[i * (col + 1) + col];
        full_cpy[i * (col + 1) + col] = full[i * (col + 1) + col];
    }

    int r1 = Rang(full, row, col + 1);
    int r2 = Rang(mat, row, col);



    if (r1 != r2)
    {
        std::cout << "NO" << std::endl;
    }
    else if (r1 < col)
    {
        std::cout << "INF" << std::endl;
    }
    else
    {
        std::cout << "YES" << std::endl;

        Gauss(full_cpy, row, col);
        for (int i = 0; i < col; i++)
        {
            std::cout << full_cpy[i * (col + 1) + col] / full_cpy[i * (col + 1) + i] << " ";
        }
    }


    std::cin.get();

    return 0;
}