#include <sgpp/base/operation/hash/common/basis/PsiHermiteBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/ZetaHermiteBasis.hpp>


#include <iostream>
#include <string>

using namespace std;
int main() {
  std::unique_ptr<sgpp::base::PsiHermiteBasis<unsigned int, unsigned int>> psibasis;
  psibasis = std::unique_ptr<sgpp::base::PsiHermiteBasis<unsigned int, unsigned int>>(
      new sgpp::base::PsiHermiteBasis<unsigned int, unsigned int>());

  std::unique_ptr<sgpp::base::ZetaHermiteBasis<unsigned int, unsigned int>> zetabasis;
  zetabasis = std::unique_ptr<sgpp::base::ZetaHermiteBasis<unsigned int, unsigned int>>(
      new sgpp::base::ZetaHermiteBasis<unsigned int, unsigned int>());

  cout << psibasis->eval(2, 3, 0.75) << "  " << zetabasis->eval(2, 3, 0.75) << endl;
  cout << psibasis->eval(2, 1, 0.25) << "  " << zetabasis->eval(2, 1, 0.25) << endl;
  cout << psibasis->eval(2, 3, 1) << "  " << zetabasis->eval(2, 3, 1) << endl;
  cout << psibasis->eval(2, 3, 0.85) << "  " << zetabasis->eval(2, 3, 0.85) << endl;

}
