#pragma once
#ifndef MD_MESH_DECEMATOR_H
#define MD_MESH_DECEMATOR_H
#include <queue>
#include <set>
#include <vector>
#include <limits>
#include<math.h>
#include<iostream>
#include<stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sstream>
#include <algorithm>

#define SARRAY_DEFAULT_MIN_SIZE 16

namespace MeshDecimation
{
    template < typename T, size_t N > class SArray
    {
    public:
        T& operator[](size_t i)
        {
            T* const data = Data();
            return data[i];
        }
        const T& operator[](size_t i) const
        {
            const T* const data = Data();
            return data[i];
        }
        size_t                  Size() const
        {
            return m_size;
        }
        T* const               Data()
        {
            return (m_maxSize == N) ? m_data0 : m_data;
        }
        const T* const         Data() const
        {
            return (m_maxSize == N) ? m_data0 : m_data;
        }
        void                    Clear()
        {
            m_size = 0;
            delete[] m_data;
            m_data = 0;
            m_maxSize = N;
        }
        void                    PopBack()
        {
            --m_size;
        }
        void                    Resize(size_t size)
        {
            if (size > m_maxSize)
            {
                T* temp = new T[size];
                memcpy(temp, Data(), m_size * sizeof(T));
                delete[] m_data;
                m_data = temp;
                m_maxSize = size;
            }
        }
        void                    PushBack(const T& value)
        {
            if (m_size == m_maxSize)
            {
                size_t maxSize = (m_maxSize << 1);
                T* temp = new T[maxSize];
                memcpy(temp, Data(), m_maxSize * sizeof(T));
                delete[] m_data;
                m_data = temp;
                m_maxSize = maxSize;
            }
            T* const data = Data();
            data[m_size++] = value;
        }
        bool                    Find(const T& value, size_t& pos)
        {
            T* const data = Data();
            for (pos = 0; pos < m_size; ++pos)
                if (value == data[pos]) return true;
            return false;
        }
        bool                    Insert(const T& value)
        {
            size_t pos;
            if (Find(value, pos)) return false;
            PushBack(value);
            return true;
        }
        bool                    Erase(const T& value)
        {
            size_t pos;
            T* const data = Data();
            if (Find(value, pos))
            {
                for (size_t j = pos + 1; j < m_size; ++j)
                    data[j - 1] = data[j];
                --m_size;
                return true;
            }
            return false;
        }
        void                    operator=(const SArray& rhs)
        {
            if (m_maxSize < rhs.m_size)
            {
                delete[] m_data;
                m_maxSize = rhs.m_maxSize;
                m_data = new T[m_maxSize];
            }
            m_size = rhs.m_size;
            memcpy(Data(), rhs.Data(), m_size * sizeof(T));
        }
        void                    Initialize()
        {
            m_data = 0;
            m_size = 0;
            m_maxSize = N;
        }
        SArray(const SArray& rhs)
        {
            m_data = 0;
            m_size = 0;
            m_maxSize = N;
            *this = rhs;
        }
        SArray()
        {
            Initialize();
        }
        ~SArray()
        {
            delete[] m_data;
        }
    private:
        T                       m_data0[N];
        T* m_data;
        size_t                  m_size;
        size_t                  m_maxSize;
    };

    template < typename T > class Vec3
    {
    public:
        T& operator[](size_t i) { return m_data[i]; }
        const T& operator[](size_t i) const { return m_data[i]; }
        T& X();
        T& Y();
        T& Z();
        const T& X() const;
        const T& Y() const;
        const T& Z() const;
        void				Normalize();
        T					GetNorm() const;
        void				operator= (const Vec3& rhs);
        void				operator+=(const Vec3& rhs);
        void				operator-=(const Vec3& rhs);
        void				operator-=(T a);
        void				operator+=(T a);
        void				operator/=(T a);
        void				operator*=(T a);
        Vec3				operator^ (const Vec3& rhs) const;
        T			    	operator* (const Vec3& rhs) const;
        Vec3				operator+ (const Vec3& rhs) const;
        Vec3				operator- (const Vec3& rhs) const;
        Vec3				operator- () const;
        Vec3				operator* (T rhs) const;
        Vec3				operator/ (T rhs) const;
        Vec3();
        Vec3(T a);
        Vec3(T x, T y, T z);
        Vec3(const Vec3& rhs);
        /*virtual*/			~Vec3(void);

    private:
        T					m_data[3];
    };
    template<typename T>
    const bool Colinear(const Vec3<T>& a, const Vec3<T>& b, const Vec3<T>& c);
    template<typename T>
    const T Volume(const Vec3<T>& a, const Vec3<T>& b, const Vec3<T>& c, const Vec3<T>& d);

    template <typename T>
    inline Vec3<T> operator*(T lhs, const Vec3<T>& rhs)
    {
        return Vec3<T>(lhs * rhs.X(), lhs * rhs.Y(), lhs * rhs.Z());
    }
    template <typename T>
    inline T& Vec3<T>::X()
    {
        return m_data[0];
    }
    template <typename T>
    inline  T& Vec3<T>::Y()
    {
        return m_data[1];
    }
    template <typename T>
    inline  T& Vec3<T>::Z()
    {
        return m_data[2];
    }
    template <typename T>
    inline  const T& Vec3<T>::X() const
    {
        return m_data[0];
    }
    template <typename T>
    inline  const T& Vec3<T>::Y() const
    {
        return m_data[1];
    }
    template <typename T>
    inline  const T& Vec3<T>::Z() const
    {
        return m_data[2];
    }
    template <typename T>
    inline  void Vec3<T>::Normalize()
    {
        T n = sqrt(m_data[0] * m_data[0] + m_data[1] * m_data[1] + m_data[2] * m_data[2]);
        if (n != 0.0) (*this) /= n;
    }
    template <typename T>
    inline  T Vec3<T>::GetNorm() const
    {
        return sqrt(m_data[0] * m_data[0] + m_data[1] * m_data[1] + m_data[2] * m_data[2]);
    }
    template <typename T>
    inline  void Vec3<T>::operator= (const Vec3& rhs)
    {
        this->m_data[0] = rhs.m_data[0];
        this->m_data[1] = rhs.m_data[1];
        this->m_data[2] = rhs.m_data[2];
    }
    template <typename T>
    inline  void Vec3<T>::operator+=(const Vec3& rhs)
    {
        this->m_data[0] += rhs.m_data[0];
        this->m_data[1] += rhs.m_data[1];
        this->m_data[2] += rhs.m_data[2];
    }
    template <typename T>
    inline void Vec3<T>::operator-=(const Vec3& rhs)
    {
        this->m_data[0] -= rhs.m_data[0];
        this->m_data[1] -= rhs.m_data[1];
        this->m_data[2] -= rhs.m_data[2];
    }
    template <typename T>
    inline void Vec3<T>::operator-=(T a)
    {
        this->m_data[0] -= a;
        this->m_data[1] -= a;
        this->m_data[2] -= a;
    }
    template <typename T>
    inline void Vec3<T>::operator+=(T a)
    {
        this->m_data[0] += a;
        this->m_data[1] += a;
        this->m_data[2] += a;
    }
    template <typename T>
    inline void Vec3<T>::operator/=(T a)
    {
        this->m_data[0] /= a;
        this->m_data[1] /= a;
        this->m_data[2] /= a;
    }
    template <typename T>
    inline void Vec3<T>::operator*=(T a)
    {
        this->m_data[0] *= a;
        this->m_data[1] *= a;
        this->m_data[2] *= a;
    }
    template <typename T>
    inline Vec3<T> Vec3<T>::operator^ (const Vec3<T>& rhs) const
    {
        return Vec3<T>(m_data[1] * rhs.m_data[2] - m_data[2] * rhs.m_data[1],
            m_data[2] * rhs.m_data[0] - m_data[0] * rhs.m_data[2],
            m_data[0] * rhs.m_data[1] - m_data[1] * rhs.m_data[0]);
    }
    template <typename T>
    inline T Vec3<T>::operator*(const Vec3<T>& rhs) const
    {
        return (m_data[0] * rhs.m_data[0] + m_data[1] * rhs.m_data[1] + m_data[2] * rhs.m_data[2]);
    }
    template <typename T>
    inline Vec3<T> Vec3<T>::operator+(const Vec3<T>& rhs) const
    {
        return Vec3<T>(m_data[0] + rhs.m_data[0], m_data[1] + rhs.m_data[1], m_data[2] + rhs.m_data[2]);
    }
    template <typename T>
    inline  Vec3<T> Vec3<T>::operator-(const Vec3<T>& rhs) const
    {
        return Vec3<T>(m_data[0] - rhs.m_data[0], m_data[1] - rhs.m_data[1], m_data[2] - rhs.m_data[2]);
    }
    template <typename T>
    inline  Vec3<T> Vec3<T>::operator-() const
    {
        return Vec3<T>(-m_data[0], -m_data[1], -m_data[2]);
    }

    template <typename T>
    inline Vec3<T> Vec3<T>::operator*(T rhs) const
    {
        return Vec3<T>(rhs * this->m_data[0], rhs * this->m_data[1], rhs * this->m_data[2]);
    }
    template <typename T>
    inline Vec3<T> Vec3<T>::operator/ (T rhs) const
    {
        return Vec3<T>(m_data[0] / rhs, m_data[1] / rhs, m_data[2] / rhs);
    }
    template <typename T>
    inline Vec3<T>::Vec3(T a)
    {
        m_data[0] = m_data[1] = m_data[2] = a;
    }
    template <typename T>
    inline Vec3<T>::Vec3(T x, T y, T z)
    {
        m_data[0] = x;
        m_data[1] = y;
        m_data[2] = z;
    }
    template <typename T>
    inline Vec3<T>::Vec3(const Vec3& rhs)
    {
        m_data[0] = rhs.m_data[0];
        m_data[1] = rhs.m_data[1];
        m_data[2] = rhs.m_data[2];
    }
    template <typename T>
    inline Vec3<T>::~Vec3(void) {};

    template <typename T>
    inline Vec3<T>::Vec3() {}

    template<typename T>
    inline const bool Colinear(const Vec3<T>& a, const Vec3<T>& b, const Vec3<T>& c)
    {
        return  ((c.Z() - a.Z()) * (b.Y() - a.Y()) - (b.Z() - a.Z()) * (c.Y() - a.Y()) == 0.0 /*EPS*/) &&
            ((b.Z() - a.Z()) * (c.X() - a.X()) - (b.X() - a.X()) * (c.Z() - a.Z()) == 0.0 /*EPS*/) &&
            ((b.X() - a.X()) * (c.Y() - a.Y()) - (b.Y() - a.Y()) * (c.X() - a.X()) == 0.0 /*EPS*/);
    }

    template<typename T>
    inline const T Volume(const Vec3<T>& a, const Vec3<T>& b, const Vec3<T>& c, const Vec3<T>& d)
    {
        return (a - d) * ((b - d) ^ (c - d));
    }


    //    typedef double Float;
    typedef float Float;
    struct MDVertex
    {
        SArray<int, SARRAY_DEFAULT_MIN_SIZE>    m_edges;
        SArray<int, SARRAY_DEFAULT_MIN_SIZE>    m_triangles;
        Float                                   m_Q[10];
        // 0 1 2 3
        //   4 5 6
        //     7 8
        //       9
        bool                                    m_tag;
        bool                                    m_onBoundary;
    };

    struct MDEdge
    {
        int                                     m_v1;
        int                                     m_v2;
        double                                  m_qem;
        Vec3<Float>                             m_pos;
        bool                                    m_onBoundary;
        bool                                    m_tag;
    };
    struct MDEdgePriorityQueue
    {
        int                                     m_name;
        double                                  m_qem;
        inline    friend bool                   operator<(const MDEdgePriorityQueue& lhs, const MDEdgePriorityQueue& rhs) { return (lhs.m_qem > rhs.m_qem); }
        inline    friend bool                   operator>(const MDEdgePriorityQueue& lhs, const MDEdgePriorityQueue& rhs) { return (lhs.m_qem < rhs.m_qem); }
    };
    typedef void (*CallBackFunction)(const char* msg);

    class MeshDecimator
    {
    public:
        //! Sets the call-back function
        //! @param callBack pointer to the call-back function
        void                                    SetCallBack(CallBackFunction  callBack) { m_callBack = callBack; }
        //! Gives the call-back function
        //! @return pointer to the call-back function
        const CallBackFunction                  GetCallBack() const { return m_callBack; }

        inline void                             SetEColManifoldConstraint(bool ecolManifoldConstraint) { m_ecolManifoldConstraint = ecolManifoldConstraint; }
        inline size_t                           GetNVertices()const { return m_nVertices; };
        inline size_t                           GetNTriangles() const { return m_nTriangles; };
        inline size_t                           GetNEdges() const { return m_nEdges; };
        void                                    GetMeshData(Vec3<Float>* points, Vec3<int>* triangles) const;
        void                                    ReleaseMemory();
        void                                    Initialize(size_t nVertices, size_t nTriangles,
            Vec3<Float>* points,
            Vec3<int>* triangles);
        bool                                    Decimate(size_t targetNVertices = 100,
            size_t targetNTriangles = 0,
            double targetError = std::numeric_limits<double>::max());

        MeshDecimator(void);
        ~MeshDecimator(void);
    private:
        void                                    EdgeCollapse(int v1, int v2);
        int                                     GetTriangle(int v1, int v2, int v3) const;
        int                                     GetEdge(int v1, int v2) const;
        int                                     IsBoundaryEdge(int v1, int v2) const;
        bool                                    IsBoundaryVertex(int v) const;
        void                                    InitializePriorityQueue();
        void                                    InitializeQEM();
        bool                                    ManifoldConstraint(int v1, int v2) const;
        double                                  ComputeEdgeCost(int v1, int v2, Vec3<Float>& pos) const;
        bool                                    EdgeCollapse(double& error);
    private:
        Vec3<int>* m_triangles;
        Vec3<Float>* m_points;
        size_t                                  m_nPoints;
        size_t                                  m_nInitialTriangles;
        size_t                                  m_nVertices;
        size_t                                  m_nTriangles;
        size_t                                  m_nEdges;
        double                                  m_diagBB;
        std::vector<MDVertex>                   m_vertices;
        std::vector<MDEdge>                     m_edges;
        std::priority_queue<
            MDEdgePriorityQueue,
            std::vector<MDEdgePriorityQueue>,
            std::less<MDEdgePriorityQueue> >   m_pqueue;
        CallBackFunction                        m_callBack;                    //>! call-back function
        bool* m_trianglesTags;
        bool                                    m_ecolManifoldConstraint;
    };

    MeshDecimator::MeshDecimator(void)
    {
        m_triangles = 0;
        m_points = 0;
        m_nPoints = 0;
        m_nInitialTriangles = 0;
        m_nVertices = 0;
        m_nTriangles = 0;
        m_nEdges = 0;
        m_trianglesTags = 0;
        m_ecolManifoldConstraint = true;
        m_callBack = 0;
    }

    MeshDecimator::~MeshDecimator(void)
    {
        ReleaseMemory();
    }
    void MeshDecimator::ReleaseMemory()
    {
        delete[] m_trianglesTags;
        std::vector< MDVertex > emptyVertices(0);
        m_vertices.swap(emptyVertices);
        std::vector<MDEdge> emptyEdges(0);
        m_edges.swap(emptyEdges);
        m_pqueue = std::priority_queue<
            MDEdgePriorityQueue,
            std::vector<MDEdgePriorityQueue>,
            std::less<MDEdgePriorityQueue> >();
        m_triangles = 0;
        m_points = 0;
        m_nPoints = 0;
        m_nInitialTriangles = 0;
        m_nVertices = 0;
        m_nTriangles = 0;
        m_nEdges = 0;
        m_trianglesTags = 0;
    }
    void MeshDecimator::Initialize(size_t nVertices, size_t nTriangles, Vec3<Float>* points, Vec3<int>* triangles)
    {
        m_nVertices = nVertices;
        m_nTriangles = nTriangles;
        m_nInitialTriangles = nTriangles;
        m_points = points;
        m_nPoints = nVertices;
        m_triangles = triangles;
        m_trianglesTags = new bool[m_nTriangles];
        m_edges.reserve(3 * m_nTriangles);
        m_vertices.resize(m_nVertices);
        for (size_t v = 0; v < m_nVertices; ++v)
        {
            m_vertices[v].m_tag = true;
        }
        int tri[3];
        MDEdge edge;
        edge.m_tag = true;
        edge.m_onBoundary = true;
        int nEdges = 0;
        int idEdge;
        int nTris = static_cast<int>(m_nTriangles);
        for (int t = 0; t < nTris; ++t)
        {
            tri[0] = m_triangles[t].X();
            tri[1] = m_triangles[t].Y();
            tri[2] = m_triangles[t].Z();
            m_trianglesTags[t] = true;
            for (int k = 0; k < 3; ++k)
            {
                edge.m_v1 = tri[k];
                edge.m_v2 = tri[(k + 1) % 3];
                m_vertices[edge.m_v1].m_triangles.Insert(t);
                idEdge = GetEdge(edge.m_v1, edge.m_v2);
                if (idEdge == -1)
                {
                    m_edges.push_back(edge);
                    m_vertices[edge.m_v1].m_edges.Insert(nEdges);
                    m_vertices[edge.m_v2].m_edges.Insert(nEdges);
                    ++nEdges;
                }
                else
                {
                    m_edges[idEdge].m_onBoundary = false;
                }
            }
        }
        m_nEdges = static_cast<size_t>(nEdges);
        for (size_t v = 0; v < m_nVertices; ++v)
        {
            m_vertices[v].m_onBoundary = false;
            for (size_t itE = 0; itE < m_vertices[v].m_edges.Size(); ++itE)
            {
                idEdge = m_vertices[v].m_edges[itE];
                if (m_edges[idEdge].m_onBoundary)
                {
                    m_vertices[v].m_onBoundary = true;
                    break;
                }
            }
        }
    }
    int MeshDecimator::GetTriangle(int v1, int v2, int v3) const
    {
        int i, j, k;
        int idTriangle;
        for (size_t it = 0; it < m_vertices[v1].m_triangles.Size(); ++it)
        {
            idTriangle = m_vertices[v1].m_triangles[it];
            i = m_triangles[idTriangle].X();
            j = m_triangles[idTriangle].Y();
            k = m_triangles[idTriangle].Z();
            if ((i == v1 && j == v2 && k == v3) || (i == v1 && j == v3 && k == v2) ||
                (i == v2 && j == v1 && k == v3) || (i == v2 && j == v3 && k == v1) ||
                (i == v3 && j == v2 && k == v1) || (i == v3 && j == v1 && k == v2))
            {
                return idTriangle;
            }
        }
        return -1;
    }
    int MeshDecimator::GetEdge(int v1, int v2) const
    {
        int idEdge;
        for (size_t it = 0; it < m_vertices[v1].m_edges.Size(); ++it)
        {
            idEdge = m_vertices[v1].m_edges[it];
            if ((m_edges[idEdge].m_v1 == v1 && m_edges[idEdge].m_v2 == v2) ||
                (m_edges[idEdge].m_v1 == v2 && m_edges[idEdge].m_v2 == v1))
            {
                return idEdge;
            }
        }
        return -1;
    }
    void MeshDecimator::EdgeCollapse(int v1, int v2)
    {
        int u, w;
        int shift;
        int idTriangle;
        for (size_t itT = 0; itT < m_vertices[v2].m_triangles.Size(); ++itT)
        {
            idTriangle = m_vertices[v2].m_triangles[itT];
            if (m_triangles[idTriangle].X() == v2)
            {
                shift = 0;
                u = m_triangles[idTriangle].Y();
                w = m_triangles[idTriangle].Z();
            }
            else if (m_triangles[idTriangle].Y() == v2)
            {
                shift = 1;
                u = m_triangles[idTriangle].X();
                w = m_triangles[idTriangle].Z();
            }
            else
            {
                shift = 2;
                u = m_triangles[idTriangle].X();
                w = m_triangles[idTriangle].Y();
            }

            if ((u == v1) || (w == v1))
            {
                m_trianglesTags[idTriangle] = false;
                m_vertices[u].m_triangles.Erase(idTriangle);
                m_vertices[w].m_triangles.Erase(idTriangle);
                m_nTriangles--;
            }
            else if (GetTriangle(v1, u, w) == -1)
            {
                m_vertices[v1].m_triangles.Insert(idTriangle);
                m_triangles[idTriangle][shift] = v1;
            }
            else
            {
                m_trianglesTags[idTriangle] = false;
                m_vertices[u].m_triangles.Erase(idTriangle);
                m_vertices[w].m_triangles.Erase(idTriangle);
                m_nTriangles--;
            }
        }
        int idEdge;
        for (size_t itE = 0; itE < m_vertices[v2].m_edges.Size(); ++itE)
        {
            idEdge = m_vertices[v2].m_edges[itE];
            w = (m_edges[idEdge].m_v1 == v2) ? m_edges[idEdge].m_v2 : m_edges[idEdge].m_v1;
            if (w == v1)
            {
                m_edges[idEdge].m_tag = false;
                m_vertices[w].m_edges.Erase(idEdge);
                m_nEdges--;
            }
            else if (GetEdge(v1, w) == -1)
            {
                if (m_edges[idEdge].m_v1 == v2)    m_edges[idEdge].m_v1 = v1;
                else                            m_edges[idEdge].m_v2 = v1;
                m_vertices[v1].m_edges.Insert(idEdge);
            }
            else
            {
                m_edges[idEdge].m_tag = false;
                m_vertices[w].m_edges.Erase(idEdge);
                m_nEdges--;
            }
        }
        m_vertices[v2].m_tag = false;
        m_nVertices--;
        // update boundary edges
        SArray<int, 64> incidentVertices;
        incidentVertices.PushBack(v1);
        for (size_t itE = 0; itE < m_vertices[v1].m_edges.Size(); ++itE)
        {
            incidentVertices.PushBack((m_edges[idEdge].m_v1 != v1) ? m_edges[idEdge].m_v1 : m_edges[idEdge].m_v2);
            idEdge = m_vertices[v1].m_edges[itE];
            m_edges[idEdge].m_onBoundary = (IsBoundaryEdge(m_edges[idEdge].m_v1, m_edges[idEdge].m_v2) != -1);
        }
        // update boundary vertices
        int idVertex;
        for (size_t itV = 0; itV < incidentVertices.Size(); ++itV)
        {
            idVertex = incidentVertices[itV];
            m_vertices[idVertex].m_onBoundary = false;
            for (size_t itE = 0; itE < m_vertices[idVertex].m_edges.Size(); ++itE)
            {
                idEdge = m_vertices[idVertex].m_edges[itE];
                if (m_edges[idEdge].m_onBoundary)
                {
                    m_vertices[idVertex].m_onBoundary = true;
                    break;
                }
            }
        }
    }
    int MeshDecimator::IsBoundaryEdge(int v1, int v2) const
    {
        int commonTri = -1;
        int itTriangle1, itTriangle2;
        for (size_t itT1 = 0; itT1 < m_vertices[v1].m_triangles.Size(); ++itT1)
        {
            itTriangle1 = m_vertices[v1].m_triangles[itT1];
            for (size_t itT2 = 0; itT2 < m_vertices[v2].m_triangles.Size(); ++itT2)
            {
                itTriangle2 = m_vertices[v2].m_triangles[itT2];
                if (itTriangle1 == itTriangle2)
                {
                    if (commonTri == -1)
                    {
                        commonTri = itTriangle1;
                    }
                    else
                    {
                        return -1;
                    }
                }
            }
        }
        return commonTri;
    }
    bool MeshDecimator::IsBoundaryVertex(int v) const
    {
        int idEdge;
        for (size_t itE = 0; itE < m_vertices[v].m_edges.Size(); ++itE)
        {
            idEdge = m_vertices[v].m_edges[itE];
            if (IsBoundaryEdge(m_edges[idEdge].m_v1, m_edges[idEdge].m_v2) != -1) return true;
        }
        return false;
    }

    void MeshDecimator::GetMeshData(Vec3<Float>* points, Vec3<int>* triangles) const
    {
        int* map = new int[m_nPoints];
        int counter = 0;
        for (size_t v = 0; v < m_nPoints; ++v)
        {
            if (m_vertices[v].m_tag)
            {
                points[counter] = m_points[v];
                map[v] = counter++;
            }
        }
        counter = 0;
        for (size_t t = 0; t < m_nInitialTriangles; ++t)
        {
            if (m_trianglesTags[t])
            {
                triangles[counter].X() = map[m_triangles[t].X()];
                triangles[counter].Y() = map[m_triangles[t].Y()];
                triangles[counter].Z() = map[m_triangles[t].Z()];
                counter++;
            }
        }
        delete[] map;
    }

    void MeshDecimator::InitializeQEM()
    {
        Vec3<Float> coordMin = m_points[0];
        Vec3<Float> coordMax = m_points[0];
        Vec3<Float> coord;
        for (size_t p = 1; p < m_nPoints; ++p)
        {
            coord = m_points[p];
            if (coordMin.X() > coord.X()) coordMin.X() = coord.X();
            if (coordMin.Y() > coord.Y()) coordMin.Y() = coord.Y();
            if (coordMin.Z() > coord.Z()) coordMin.Z() = coord.Z();
            if (coordMax.X() < coord.X()) coordMax.X() = coord.X();
            if (coordMax.Y() < coord.Y()) coordMax.Y() = coord.Y();
            if (coordMax.Z() < coord.Z()) coordMax.Z() = coord.Z();
        }
        coordMax -= coordMin;
        m_diagBB = coordMax.GetNorm();

        int i, j, k;
        Vec3<Float> n;
        Float d = 0;
        Float area = 0;
        for (size_t v = 0; v < m_nPoints; ++v)
        {
            memset(m_vertices[v].m_Q, 0, 10 * sizeof(Float));
            int idTriangle;
            for (size_t itT = 0; itT < m_vertices[v].m_triangles.Size(); ++itT)
            {
                idTriangle = m_vertices[v].m_triangles[itT];
                i = m_triangles[idTriangle].X();
                j = m_triangles[idTriangle].Y();
                k = m_triangles[idTriangle].Z();
                n = (m_points[j] - m_points[i]) ^ (m_points[k] - m_points[i]);
                area = n.GetNorm();
                n.Normalize();
                d = -(m_points[v] * n);
                m_vertices[v].m_Q[0] += area * (n.X() * n.X());
                m_vertices[v].m_Q[1] += area * (n.X() * n.Y());
                m_vertices[v].m_Q[2] += area * (n.X() * n.Z());
                m_vertices[v].m_Q[3] += area * (n.X() * d);
                m_vertices[v].m_Q[4] += area * (n.Y() * n.Y());
                m_vertices[v].m_Q[5] += area * (n.Y() * n.Z());
                m_vertices[v].m_Q[6] += area * (n.Y() * d);
                m_vertices[v].m_Q[7] += area * (n.Z() * n.Z());
                m_vertices[v].m_Q[8] += area * (n.Z() * d);
                m_vertices[v].m_Q[9] += area * (d * d);
            }
        }
        Vec3<Float> u1, u2;
        const Float w = static_cast<Float>(1000);
        int t, v1, v2, v3;
        for (size_t e = 0; e < m_edges.size(); ++e)
        {
            v1 = m_edges[e].m_v1;
            v2 = m_edges[e].m_v2;
            t = IsBoundaryEdge(v1, v2);
            if (t != -1)
            {
                if (m_triangles[t].X() != v1 && m_triangles[t].X() != v2) v3 = m_triangles[t].X();
                else if (m_triangles[t].Y() != v1 && m_triangles[t].Y() != v2) v3 = m_triangles[t].Y();
                else                                                           v3 = m_triangles[t].Z();
                u1 = m_points[v2] - m_points[v1];
                u2 = m_points[v3] - m_points[v1];
                area = w * (u1 ^ u2).GetNorm();
                u1.Normalize();
                n = u2 - (u2 * u1) * u1;
                n.Normalize();

                d = -(m_points[v1] * n);
                m_vertices[v1].m_Q[0] += area * (n.X() * n.X());
                m_vertices[v1].m_Q[1] += area * (n.X() * n.Y());
                m_vertices[v1].m_Q[2] += area * (n.X() * n.Z());
                m_vertices[v1].m_Q[3] += area * (n.X() * d);
                m_vertices[v1].m_Q[4] += area * (n.Y() * n.Y());
                m_vertices[v1].m_Q[5] += area * (n.Y() * n.Z());
                m_vertices[v1].m_Q[6] += area * (n.Y() * d);
                m_vertices[v1].m_Q[7] += area * (n.Z() * n.Z());
                m_vertices[v1].m_Q[8] += area * (n.Z() * d);
                m_vertices[v1].m_Q[9] += area * (d * d);

                d = -(m_points[v2] * n);
                m_vertices[v2].m_Q[0] += area * (n.X() * n.X());
                m_vertices[v2].m_Q[1] += area * (n.X() * n.Y());
                m_vertices[v2].m_Q[2] += area * (n.X() * n.Z());
                m_vertices[v2].m_Q[3] += area * (n.X() * d);
                m_vertices[v2].m_Q[4] += area * (n.Y() * n.Y());
                m_vertices[v2].m_Q[5] += area * (n.Y() * n.Z());
                m_vertices[v2].m_Q[6] += area * (n.Y() * d);
                m_vertices[v2].m_Q[7] += area * (n.Z() * n.Z());
                m_vertices[v2].m_Q[8] += area * (n.Z() * d);
                m_vertices[v2].m_Q[9] += area * (d * d);
            }
        }
    }
    void MeshDecimator::InitializePriorityQueue()
    {
        double progressOld = -1.0;
        double progress = 0.0;
        char msg[1024];
        double ptgStep = 1.0;
        int v1, v2;
        MDEdgePriorityQueue pqEdge;
        size_t nE = m_edges.size();
        for (size_t e = 0; e < nE; ++e)
        {
            progress = e * 100.0 / nE;
            if (fabs(progress - progressOld) > ptgStep && m_callBack)
            {
                sprintf(msg, "%3.2f %% \t \t \r", progress);
                (*m_callBack)(msg);
                progressOld = progress;
            }

            if (m_edges[e].m_tag)
            {
                v1 = m_edges[e].m_v1;
                v2 = m_edges[e].m_v2;
                if ((!m_ecolManifoldConstraint) || (ManifoldConstraint(v1, v2)))
                {
                    pqEdge.m_qem = m_edges[e].m_qem = ComputeEdgeCost(v1, v2, m_edges[e].m_pos);
                    pqEdge.m_name = static_cast<int>(e);
                    m_pqueue.push(pqEdge);
                }
            }
        }
    }
    double MeshDecimator::ComputeEdgeCost(int v1, int v2, Vec3<Float>& newPos) const
    {
        double Q[10];
        double M[12];
        Vec3<double> pos;
        for (int i = 0; i < 10; ++i) Q[i] = m_vertices[v1].m_Q[i] + m_vertices[v2].m_Q[i];
        M[0] = Q[0]; // (0, 0)
        M[1] = Q[1]; // (0, 1) 
        M[2] = Q[2]; // (0, 2)
        M[3] = Q[3]; // (0, 3)
        M[4] = Q[1]; // (1, 0)
        M[5] = Q[4]; // (1, 1)
        M[6] = Q[5]; // (1, 2)
        M[7] = Q[6]; // (1, 3)
        M[8] = Q[2]; // (2, 0)
        M[9] = Q[5]; // (2, 1)
        M[10] = Q[7]; // (2, 2);
        M[11] = Q[8]; // (2, 3);
        double det = M[0] * M[5] * M[10] + M[1] * M[6] * M[8] + M[2] * M[4] * M[9]
            - M[0] * M[6] * M[9] - M[1] * M[4] * M[10] - M[2] * M[5] * M[8];
        if (det != 0.0)
        {
            double d = 1.0 / det;
            pos.X() = d * (M[1] * M[7] * M[10] + M[2] * M[5] * M[11] + M[3] * M[6] * M[9]
                - M[1] * M[6] * M[11] - M[2] * M[7] * M[9] - M[3] * M[5] * M[10]);
            pos.Y() = d * (M[0] * M[6] * M[11] + M[2] * M[7] * M[8] + M[3] * M[4] * M[10]
                - M[0] * M[7] * M[10] - M[2] * M[4] * M[11] - M[3] * M[6] * M[8]);
            pos.Z() = d * (M[0] * M[7] * M[9] + M[1] * M[4] * M[11] + M[3] * M[5] * M[8]
                - M[0] * M[5] * M[11] - M[1] * M[7] * M[8] - M[3] * M[4] * M[9]);
            newPos.X() = static_cast<Float>(pos.X());
            newPos.Y() = static_cast<Float>(pos.Y());
            newPos.Z() = static_cast<Float>(pos.Z());
        }
        else
        {
            const Float w = static_cast<Float>(0.5f);
            newPos = w * m_points[v1] + w * m_points[v2];
            pos.X() = static_cast<double>(newPos.X());
            pos.Y() = static_cast<double>(newPos.Y());
            pos.Z() = static_cast<double>(newPos.Z());
        }

        double qem = pos.X() * (Q[0] * pos.X() + Q[1] * pos.Y() + Q[2] * pos.Z() + Q[3]) +
            pos.Y() * (Q[1] * pos.X() + Q[4] * pos.Y() + Q[5] * pos.Z() + Q[6]) +
            pos.Z() * (Q[2] * pos.X() + Q[5] * pos.Y() + Q[7] * pos.Z() + Q[8]) +
            (Q[3] * pos.X() + Q[6] * pos.Y() + Q[8] * pos.Z() + Q[9]);

        Vec3<Float> d1;
        Vec3<Float> d2;
        Vec3<Float> n1;
        Vec3<Float> n2;
        Vec3<Float> oldPosV1 = m_points[v1];
        Vec3<Float> oldPosV2 = m_points[v2];

        SArray<int, SARRAY_DEFAULT_MIN_SIZE> triangles = m_vertices[v1].m_triangles;
        int idTriangle;
        for (size_t itT = 0; itT < m_vertices[v2].m_triangles.Size(); ++itT)
        {
            idTriangle = m_vertices[v2].m_triangles[itT];
            triangles.Insert(idTriangle);
        }
        int a[3];
        for (size_t itT = 0; itT != triangles.Size(); ++itT)
        {
            idTriangle = triangles[itT];
            a[0] = m_triangles[idTriangle].X();
            a[1] = m_triangles[idTriangle].Y();
            a[2] = m_triangles[idTriangle].Z();

            d1 = m_points[a[1]] - m_points[a[0]];
            d2 = m_points[a[2]] - m_points[a[0]];
            n1 = d1 ^ d2;

            m_points[v1] = newPos;
            m_points[v2] = newPos;

            d1 = m_points[a[1]] - m_points[a[0]];
            d2 = m_points[a[2]] - m_points[a[0]];
            n2 = d1 ^ d2;

            m_points[v1] = oldPosV1;
            m_points[v2] = oldPosV2;

            n1.Normalize();
            n2.Normalize();
            if (n1 * n2 < 0.0)
            {
                return std::numeric_limits<double>::max();
            }
        }
        if (m_ecolManifoldConstraint && !ManifoldConstraint(v1, v2))
        {
            return std::numeric_limits<double>::max();
        }
        return qem;
    }
    bool MeshDecimator::ManifoldConstraint(int v1, int v2) const
    {
        std::set<int> vertices;
        int a, b;
        int idEdge1;
        int idEdge2;
        int idEdgeV1V2;
        for (size_t itE1 = 0; itE1 < m_vertices[v1].m_edges.Size(); ++itE1)
        {
            idEdge1 = m_vertices[v1].m_edges[itE1];
            a = (m_edges[idEdge1].m_v1 == v1) ? m_edges[idEdge1].m_v2 : m_edges[idEdge1].m_v1;
            vertices.insert(a);
            if (a != v2)
            {
                for (size_t itE2 = 0; itE2 < m_vertices[v2].m_edges.Size(); ++itE2)
                {
                    idEdge2 = m_vertices[v2].m_edges[itE2];
                    b = (m_edges[idEdge2].m_v1 == v2) ? m_edges[idEdge2].m_v2 : m_edges[idEdge2].m_v1;
                    vertices.insert(b);
                    if (a == b)
                    {
                        if (GetTriangle(v1, v2, a) == -1)
                        {
                            return false;
                        }
                    }
                }
            }
            else
            {
                idEdgeV1V2 = idEdge1;
            }
        }
        if (vertices.size() <= 4 || (m_vertices[v1].m_onBoundary && m_vertices[v2].m_onBoundary && !m_edges[idEdgeV1V2].m_onBoundary))
        {
            return false;
        }
        return true;
    }
    bool MeshDecimator::EdgeCollapse(double& qem)
    {
        MDEdgePriorityQueue currentEdge;
        int v1, v2;
        bool done = false;
        do
        {
            done = false;
            if (m_pqueue.size() == 0)
            {
                done = true;
                break;
            }
            else
            {
                currentEdge = m_pqueue.top();
                m_pqueue.pop();
            }
        } while ((!m_edges[currentEdge.m_name].m_tag) || (m_edges[currentEdge.m_name].m_qem != currentEdge.m_qem));

        if (done) return false;
        v1 = m_edges[currentEdge.m_name].m_v1;
        v2 = m_edges[currentEdge.m_name].m_v2;

        qem = currentEdge.m_qem;
        EdgeCollapse(v1, v2);
        m_points[v1] = m_edges[currentEdge.m_name].m_pos;
        for (int k = 0; k < 10; k++) m_vertices[v1].m_Q[k] += m_vertices[v2].m_Q[k];

        // Update priority queue
        int idEdge;
        int a, b;
        SArray<int, SARRAY_DEFAULT_MIN_SIZE> incidentVertices;
        for (size_t itE = 0; itE < m_vertices[v1].m_edges.Size(); ++itE)
        {
            idEdge = m_vertices[v1].m_edges[itE];
            a = m_edges[idEdge].m_v1;
            b = m_edges[idEdge].m_v2;
            incidentVertices.PushBack((a != v1) ? a : b);
            MDEdgePriorityQueue pqEdge;
            pqEdge.m_qem = m_edges[idEdge].m_qem = ComputeEdgeCost(a, b, m_edges[idEdge].m_pos);
            pqEdge.m_name = idEdge;
            m_pqueue.push(pqEdge);
        }
        int idVertex;
        for (size_t itV = 0; itV < incidentVertices.Size(); ++itV)
        {
            idVertex = incidentVertices[itV];
            for (size_t itE = 0; itE < m_vertices[idVertex].m_edges.Size(); ++itE)
            {
                idEdge = m_vertices[idVertex].m_edges[itE];
                a = m_edges[idEdge].m_v1;
                b = m_edges[idEdge].m_v2;
                if (a != v1 && b != v1)
                {
                    MDEdgePriorityQueue pqEdge;
                    pqEdge.m_qem = m_edges[idEdge].m_qem = ComputeEdgeCost(a, b, m_edges[idEdge].m_pos);
                    pqEdge.m_name = idEdge;
                    m_pqueue.push(pqEdge);
                }
            }
        }
        return true;
    }
    bool MeshDecimator::Decimate(size_t targetNVertices, size_t targetNTriangles, double targetError)
    {
        double qem = 0.0;
        double progressOld = -1.0;
        double progress = 0.0;
        char msg[1024];
        double ptgStep = 1.0;

        if (m_callBack)
        {
            std::ostringstream msg;
            msg << "+ Mesh" << std::endl;
            msg << "\t # vertices                     \t" << m_nPoints << std::endl;
            msg << "\t # triangles                    \t" << m_nTriangles << std::endl;
            msg << "+ Parameters" << std::endl;
            msg << "\t target # of vertices           \t" << targetNVertices << std::endl;
            msg << "\t target # of triangles          \t" << targetNTriangles << std::endl;
            msg << "\t QEM                            \t" << targetError << std::endl;
            (*m_callBack)(msg.str().c_str());
        }

        if (m_callBack) (*m_callBack)("+ Initialize QEM \n");
        InitializeQEM();
        if (m_callBack) (*m_callBack)("+ Initialize priority queue \n");
        InitializePriorityQueue();
        if (m_callBack) (*m_callBack)("+ Simplification \n");
        double invDiag = 1.0 / m_diagBB;
        while ((m_pqueue.size() > 0) &&
            (m_nEdges > 0) &&
            (m_nVertices > targetNVertices) &&
            (m_nTriangles > targetNTriangles) &&
            (qem < targetError))
        {
            progress = 100.0 - m_nVertices * 100.0 / m_nPoints;
            if (fabs(progress - progressOld) > ptgStep && m_callBack)
            {
                sprintf(msg, "%3.2f %% V = %lu \t QEM = %f \t \t \r", progress, static_cast<unsigned int>(m_nVertices), qem);
                (*m_callBack)(msg);
                progressOld = progress;
            }
            if (!EdgeCollapse(qem)) break;
            if (qem < 0.0) qem = 0.0;
            else           qem = sqrt(qem) * invDiag;
        }
        if (m_callBack)
        {
            std::ostringstream msg;
            msg << "+ Simplification output" << std::endl;
            msg << "\t # vertices                     \t" << m_nVertices << std::endl;
            msg << "\t # triangles                    \t" << m_nTriangles << std::endl;
            msg << "\t QEM                            \t" << qem << std::endl;
            (*m_callBack)(msg.str().c_str());
        }
        return true;
    }
}
#endif