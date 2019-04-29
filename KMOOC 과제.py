import ast
def transpose_matrix(m):
    return list(map(list, zip(*m))) #* unpacking후 zip으로 묶어서 transpose실행 하는 transpoes matrix 함수이다.

def sub_matrix(m,i,j):
    return [row[:j] + row[j+1:] for row in (m[:i]+m[i+1:])] #i번째행의 원소와 j번째 열의 원소를 제외한 나머지 원소들로만 구성한 sub matrix를 구하는 서브함수이다.

def Determinant_matrix(m): # 행렬식을 구하는 서브함수입니다.

    if len(m) == 2:
        return m[0][0]*m[1][1]-m[0][1]*m[1][0] #N=2일때 ad-bc인 공식을 사용한다.

    determinant = 0
    for c in range(len(m)):
        determinant += ((-1)**c)*m[0][c]*Determinant_matrix(sub_matrix(m,0,c)) # 행렬식 공식인 'Sigma_{j=1}^{n} (-1)^(j+1) a_{ij} detA_{1j}'를 코드로 작성하였습니다.
    return determinant

def minor_matrix(m,r,c): # minor_matrix를 구하는 함수
    return Determinant_matrix(sub_matrix(m,r,c))

def adjoint_matrix(m): # adjoint matrix를 구하는 서브함수입니다.
    cofactors = []
    for r in range(len(m)):
        cofactorRow = []
        for c in range(len(m)):
            cofactorRow.append(((-1)**(r+c)) * minor_matrix(m,r,c))  # r번째 행의 cofactor들을 공식인 '(-1)^{r+c} * 소행렬식'을 이용하여 구한다.
        cofactors.append(cofactorRow)
    adjoint = transpose_matrix(cofactors) # cofacor들로 구성된 행렬을 transpose 해서 adjoint matrix를 구한다.
    return adjoint
  
def Inverse_Matrix(m): # 역행렬을 구하는 서브함수 입니다.
    determinant = Determinant_matrix(m)

    if len(m) == 2:                                             # 정방행렬이 2x2일때 공식인 '1/(ad-bc) * [[d,-b],[-c,a]]'을 이용해 역행렬을 구한다.
        return [[m[1][1]/determinant, -1*m[0][1]/determinant],
                [-1*m[1][0]/determinant, m[0][0]/determinant]]

    adjoint = adjoint_matrix(m)
    inverse = adjoint

    for r in range(len(adjoint)):
        for c in range(len(adjoint)):
            inverse[r][c] = round(adjoint[r][c]/determinant,4) # 공식대로 adjoint matrix에 행렬식을 나누어서 역행렬을 구한다.
    return inverse

if __name__ == "__main__":
    m = input("역행렬을 구하고 싶은 정방행렬을 입력하시오\n")   #ex) [[2,5,6],[2,5,7],[2,7,4]]의 역행렬은 [[7.25, -5.5, -1.25], [-1.5, 1.0, 0.5], [-1.0, 1.0, -0.0]] 이다.
    m = ast.literal_eval(m)
    I = Inverse_Matrix(m)
    print("역행렬은 \n",I)

    