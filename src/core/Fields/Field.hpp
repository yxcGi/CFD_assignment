#pragma once

#include "FaceField.hpp"
#include "CellField.hpp"
#include "FieldOperators/Gradient.hpp"
#include <algorithm>
#include <fstream>



// 采用模板全特化，判断梯度类型
template<typename T>
struct GradientType;

template<>
struct GradientType<Scalar>
{
    using Type = Vector<Scalar>;
};

template<>
struct GradientType<Vector<Scalar>>
{
    using Type = Tensor<Scalar>;
};

// 输出场的文件类型
enum class WriteFileType
{
    TECPLOT,
};


template <typename Tp>
class Field
{
    using ULL = unsigned long long;
    // 定义0的界限
    static constexpr Scalar EPSILON = 1e-10;
public:
    Field() = delete;
    Field(const std::string& name, Mesh* mesh);
    Field(const Field<Tp>&) = delete;
    Field(Field<Tp>&&) noexcept = default;
    Field<Tp>& operator=(const Field<Tp>&) = delete;
    Field<Tp>& operator=(Field<Tp>&&) noexcept = default;
    ~Field() = default;

public:
    /* =========赋值======== */ // 与BaseField，接口一样
    // 统一赋值
    void setValue(const Tp& value);
    // 采用函数对象
    void setValue(const std::function<Tp(Scalar, Scalar, Scalar)>& func);


    // 获取器
    const FaceField<Tp>& getFaceField() const;
    const CellField<Tp>& getCellField() const;
    std::string getName() const;
    Mesh* getMesh() const;

    // 场是否有效，必须要setValue场才有效
    bool isValid() const;
    // 边界条件是否有效
    bool isBoundaryConditionValid() const;


    // cell场到face场的插值，边界面需要用到边界条件（边界条件需设置完整），调用后会更新当前场的梯度
    void cellToFace(interpolation::Scheme scheme = interpolation::Scheme::LINEAR);



    // 设置计算梯度方法
    void setGradientMethod(GradientMethod method);

    /* --------设置边界条件-------- */

    /* --------设置边界条件-------- */
    // a * φ + b * ∂φ/∂n = c
    void setBoundaryCondition(const std::string& name, Scalar a, Scalar b, const Tp& c);


    void writeToFile(const std::string& fileName, WriteFileType fileType = WriteFileType::TECPLOT) const;


private:
    // 私有接口
    void writToTecplot(const std::string& fileName, Mesh::Dimension dim) const;     // 输出成tecplot文件



private:
    FaceField<Tp> faceField_;        // 当前面场值
    CellField<Tp> cellField_;        // 当前单元场值
    CellField<Tp> cellField_0_;      // 上一步单元场的值
    CellField<decltype(Tp()* Vector<Scalar>())> cellGradientField_; // 上一步单元场梯度
    // CellField<typename GradientType<Tp>::Type> cellGradientField_; // 单元场梯度
    // CellField<decltype(Tp()* Vector<Scalar>())> cellGradientField_; // 单元场梯度
    std::string name_;
    GradientMethod gradientMethod_{ GradientMethod::GAUSS_GREEN };      // 梯度计算方法
    std::unordered_map<std::string, BoundaryCondition<Tp>> boundaryConditions_;
    Interpolation<Tp> interpolator_;     // 插值函数对象
};

#pragma region 函数实现

template<typename Tp>
inline Field<Tp>::Field(const std::string& name, Mesh* mesh)
    : faceField_(name, mesh)
    , cellField_(name, mesh)
    , cellField_0_(name, mesh)
    , cellGradientField_(name, mesh)
    , name_(name)
    , interpolator_(Interpolation<Tp>())
{
    const std::unordered_map<std::string, BoundaryPatch>& boundaryPatches = mesh->getBoundaryPatches();
    for (const auto& [name, patch] : boundaryPatches)
    {
        boundaryConditions_.emplace(name, BoundaryCondition<Tp>(patch));
    }
}

template<typename Tp>
inline void Field<Tp>::setValue(const Tp& value)
{
    faceField_.setValue(value);
    cellField_.setValue(value);
    cellField_0_.setValue(value);

    // 利用面值计算单元中心梯度
    cellGradientField_ = grad(*this, gradientMethod_);
}

template<typename Tp>
inline void Field<Tp>::setValue(const std::function<Tp(Scalar, Scalar, Scalar)>& func)
{
    faceField_.setValue(func);
    cellField_.setValue(func);
    cellField_0_.setValue(func);

    // 利用面值计算单元中心梯度
    cellGradientField_ = grad(*this, gradientMethod_);
}


template<typename Tp>
inline const FaceField<Tp>& Field<Tp>::getFaceField() const
{
    return faceField_;
}

template<typename Tp>
inline const CellField<Tp>& Field<Tp>::getCellField() const
{
    return cellField_;
}

template<typename Tp>
inline std::string Field<Tp>::getName() const
{
    return name_;
}

template<typename Tp>
inline Mesh* Field<Tp>::getMesh() const
{
    return faceField_.getMesh();
}

template<typename Tp>
inline bool Field<Tp>::isValid() const
{
    return faceField_.isValid() && cellField_.isValid();
}

template<typename Tp>
inline bool Field<Tp>::isBoundaryConditionValid() const
{
    return std::all_of(
        boundaryConditions_.begin(),
        boundaryConditions_.end(),
        [](const std::pair<std::string, BoundaryCondition<Tp>>& bc) {
            return (
                bc.second.isValid() ||
                bc.second.getType() == BoundaryPatch::BoundaryType::EMPTY
                );
        }
    );
}

template<typename Tp>
inline void Field<Tp>::cellToFace(interpolation::Scheme scheme)
{
    if (!isValid())
    {
        std::cerr << "Field is not valid!" << std::endl;
        throw std::runtime_error("Field is not valid!");
    }
    if (!isBoundaryConditionValid())
    {
        std::cerr << "Boundary condition is not valid!" << std::endl;
        throw std::runtime_error("Boundary condition is not valid!");
    }

    using LL = long long;
    Mesh* mesh = this->getMesh();      // 获取网格
    const std::vector<Face>& faces = mesh->getFaces();  // 面列表
    const std::vector<Cell>& cells = mesh->getCells();

    // 先遍历内部面
    std::vector<ULL> internalFaceIndexes = mesh->getInternalFaceIndexes();
    const CellField<Tp>& cellField = this->getCellField();

    for (const ULL internalFaceIndex : internalFaceIndexes)     // 遍历内部面
    {
        // 获取当前面，以及其相邻单元索引
        const Face& face = faces[internalFaceIndex];
        ULL ownerIndex = face.getOwnerIndex();
        LL neighborIndex = face.getNeighborIndex();

        // 获取当前面，以及其相邻单元中心坐标
        Vector<Scalar> faceCenter = face.getCenter();
        Vector<Scalar> ownerCenter = cells[ownerIndex].getCenter();
        Vector<Scalar> neighborCenter = cells[neighborIndex].getCenter();

        // 获取面的法向量
        Vector<Scalar> faceNormal = face.getNormal();

        // 计算面到两个面的距离
        Scalar ownerDistance = (faceCenter - ownerCenter) & faceNormal;
        Scalar neighborDistance = (faceCenter - neighborCenter) & faceNormal;

        // 计算插值权重
        Scalar alpha = ownerDistance / (ownerDistance + neighborDistance);

        // 插值
        Tp ownerValue = cellField[ownerIndex];
        Tp neighborValue = cellField[neighborIndex];

        Tp faceValue = interpolator_(ownerValue, neighborValue, scheme, alpha);

        // 设置面场
        faceField_.setValue(internalFaceIndex, faceValue);
    }

    // 再遍历边界面，需考虑边界条件
    for (const auto& [name, bc] : cellField_0_.getBoundaryConditions())
    {
        ULL nFace = bc.getNFace();
        ULL startFace = bc.getStartFace();

        if (bc.getType() ==
            BoundaryPatch::BoundaryType::EMPTY) // 不需要处理empty边界
        {
            continue;
        }


        for (ULL boundaryFaceIndex = startFace;
            boundaryFaceIndex < startFace + nFace;
            ++boundaryFaceIndex)
        {
            const Face& face = faces[boundaryFaceIndex];
            ULL ownerIndex = face.getOwnerIndex();
            const Cell& ownerCell = cells[ownerIndex];
            const Vector<Tp>& faceCenter = face.getCenter();
            const Vector<Tp>& ownerCenter = ownerCell.getCenter();
            const Vector<Scalar>& normal = face.getNormal();   // 面法向量
            Vector<Scalar> V_CB = faceCenter - ownerCenter;

            // 计算中间量, normal = E + T
            Vector<Scalar> E = V_CB / (V_CB & normal);
            Vector<Scalar> T = normal - E;
            Scalar E_magnitude = E.magnitude();
            Scalar distance_CB = ownerCenter.getDistance(faceCenter);
            Scalar a = bc.get_a();
            Scalar b = bc.get_b();
            const Tp& c = bc.get_c();
            // const auto& gradientCell =       // 计算梯度，挖坑


            // 计算c1, c2
            Scalar c1 = (b * E_magnitude) / (a * distance_CB + b * E_magnitude);
            // auto c2 = (c - b * )    // 需要梯度计算，挖坑
            const decltype(Tp() * Vector<Scalar>())& cellGradient = cellGradientField_[ownerIndex];
            Tp c2 = (c - (b * cellGradient & T)) * distance_CB / (a * distance_CB + b * E_magnitude);

            // 计算边界面值
            Tp boundaryFaceValue = c1 * cellField[ownerIndex] + c2;

            faceField_.setValue(boundaryFaceIndex, boundaryFaceValue);
        }
    }
    // 利用新的面值记录计算本时间步的梯度，用于下一时间步的边界面值计算
    cellGradientField_ = grad(*this);
}

template<typename Tp>
inline void Field<Tp>::setGradientMethod(GradientMethod method)
{
    gradientMethod_ = method;
}

template<typename Tp>
inline void Field<Tp>::setBoundaryCondition(const std::string& name, Scalar a, Scalar b, const Tp& c)
{
    if (!isValid())
    {
        throw std::runtime_error("Field is not valid!");
    }

    // 判断是否存在该边界
    using It = typename std::unordered_map<std::string, BoundaryCondition<Tp>>::iterator;
    It it = boundaryConditions_.find(name);
    if (it == boundaryConditions_.end())    // 没找到该边界
    {
        std::cerr << "Boundary condition " << name << " does not exist!" << std::endl;
        throw std::runtime_error("Boundary condition does not exist!");
    }
    it->second.setBoundaryCondition(name, a, b, c);
    it->second.setValid();      // 设置该边界条件有效
}

template<typename Tp>
inline void Field<Tp>::writeToFile(const std::string& fileName, WriteFileType fileType) const
{
    if (!isValid())
    {
        std::cerr << "Field is not valid! cannot write to file!" << std::endl;
        throw std::runtime_error("Field is not valid!");
    }

    Mesh::Dimension dim = this->getMesh()->getDimension();
    if (dim == Mesh::Dimension::TWO_D)   // 二维
    {
        std::cout << "The two-dimensional field is being output..." << std::endl;
        if (fileType == WriteFileType::TECPLOT)
        {
            writToTecplot(fileName, dim);
        }
    }
    else        // 三维
    {
        std::cout << "The three-dimensional field is begin output..." << std::endl;
        if (fileType == WriteFileType::TECPLOT)
        {
            writToTecplot(fileName, dim);
        }
    }
}

template<typename Tp>
inline void Field<Tp>::writToTecplot(const std::string& fileName, Mesh::Dimension dim) const
{
    // 私有接口，公有接口调用输出时已经判断过Field是否有效。

    Mesh* mesh = this->getMesh();

    if (dim == Mesh::Dimension::TWO_D)  // 二维
    {
        std::ofstream ofs(fileName, std::ios::trunc);
        if (!ofs.is_open())
        {
            std::cerr << "writToTecplot() Error: Unable to open file " << fileName << " for writing.\n";
            throw std::runtime_error(fileName + " file open error");
        }

        // 写入头信息
        ofs << "Title=" << "\"" << this->getName() << "_2D_Tecplot\"\n";
        if constexpr (std::is_same_v<Tp, Scalar>)       // 标量场
        {
            std::cout << "Writing scalar field to Tecplot file...\n";
            ofs << R"(VARIABLES="X","Y",")" << name_ << "\"\n";
            ofs << "ZONE T=\"" << name_ << "\",N=" << (mesh->getPointNumber()) << ",E=" << mesh->getCellNumber();   // 进入分支接着输出

            if (this->getMesh()->getMeshShape() ==
                Mesh::MeshShape::TRIANGLE)      // 二维三角形网格  标量
            {
                ofs << ",ET=triangle" << std::endl;     // 挖坑
                // 先输出所有的x，y，z
                for (const Point& point : mesh->getPoints())
                {
                    ofs << point.x() << " ";
                }
                ofs << "\n";
                for (const Point& point : mesh->getPoints())
                {
                    ofs << point.y() << " ";
                }
                ofs << "\n";
                for (const Point& point : mesh->getPoints())
                {
                    ofs << point.z() << " ";
                }
                ofs << "\n";

                // 输出每个cell的场值
                for (ULL i = 0; i < mesh->getCellNumber(); ++i)
                {
                    ofs << this->cellField_[i] << " ";
                }
                ofs << "\n";


                // 输出每个cell的点索引
                for (const Cell& cell : mesh->getCells())
                {
                    for (const ULL index : cell.getPointIndexes())
                    {
                        ofs << index << " ";
                    }
                    ofs << "\n";
                }
            }
            else if (this->getMesh()->getMeshShape() ==
                Mesh::MeshShape::QUADRILATERAL) // 二维四边形网格  标量
            {
                ofs << ",ZONETYPE=FEQUADRILATERAL\n";
                ofs << "DATAPACKING=BLOCK\n";
                ofs << "VARLOCATION=([1-2]=NODAL, [3]=CELLCENTERED)\n";

                // 输出所有点的x,y坐标
                for (const Point& point : mesh->getPoints())
                {
                    ofs << point.x() << "\n";
                }
                for (const Point& point : mesh->getPoints())
                {
                    ofs << point.y() << "\n";
                }
                // 输出所有单元的场值
                for (ULL i = 0; i < mesh->getCellNumber(); ++i)
                {
                    ofs << this->cellField_0_[i] << "\n";
                }
                // 每个二维单元由哪几个点构成
                for (const Cell& cell : mesh->getCells())
                {
                    ULL emptyFaceIndex = 0;     // 记录该单元的empty面的索引
                    for (const ULL faceId : cell.getFaceIndexes())
                    {
                        if (faceId >= mesh->getEmptyFaceIndexesPair().first &&
                            faceId < mesh->getEmptyFaceIndexesPair().second)
                        {
                            emptyFaceIndex = faceId;
                        }
                    }
                    const std::vector<ULL>& emptyFacePointIndexes = mesh->getFaces()[emptyFaceIndex].getPointIndexes();
                    for (const ULL pointId : emptyFacePointIndexes)
                    {
                        ofs << pointId + 1 << " ";      // tecplot索引从1开始
                    }
                    ofs << "\n";
                }
            }


        }
        else if constexpr (std::is_same_v<Tp, Vector<Scalar>>)    // 矢量场
        {
            std::cout << "Writing vector field to Tecplot file..." << std::endl;
            ofs << R"(VARIABLES="X","Y","U","V")" << std::endl;
            ofs << "ZONE T=\"" << name_ << "\",n=" << (mesh->getPointNumber()) / 2 << ",e=" << mesh->getCellNumber() << ",f=feblock";

            if (this->getMesh()->getMeshShape() ==
                Mesh::MeshShape::TRIANGLE)      // 二维三角形网格  矢量
            {
                ofs << ",et=triangle" << std::endl;

            }
            else if (this->getMesh()->getMeshShape() ==
                Mesh::MeshShape::QUADRILATERAL) // 二维四边形网格  矢量
            {
                ofs << ",et=quadrilateral" << std::endl;

            }


        }
    }
}