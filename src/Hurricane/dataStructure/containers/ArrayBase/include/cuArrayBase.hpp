/*!
 * \file ArrayBase.hpp
 * \brief Headers of the base of Array.
 * \author Rao Sihang
 * \version V2.0.0
 * \date 2022.05.02
 *
 * OpenHurricane: Open parts of Hurricane project (Highly Universal Rocket & Ramjet sImulation Codes for ANalysis and Evaluation)
 * \copyright Copyright (C) 2019-2024, Prof. Xu Xu's group at Beihang University.
 *
 * License
 *		This file is part of OpenHurricane
 *
 *		OpenHurricane is free software: you can redistribute it and/or modify it
 *		under the terms of the GNU General Public License as published by
 *		the Free Software Foundation, either version 3 of the License, or
 *		(at your option) any later version.
 *
 *		OpenHurricane is distributed in the hope that it will be useful, but WITHOUT
 *		ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *		FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *		for more details.
 *
 *		You should have received a copy of the GNU General Public License
 *		along with OpenHurricane.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 */

#pragma once
#include "cuAttributes.hpp"
#include "cuPreset.hpp"

#ifdef CUDA_PARALLEL
namespace OpenHurricane {
    /**
     * \brief CUDA 1D array struct.
     */
    template <typename Type> class cu1DArray {
    public:
        using value_type = Type;
        using reference = Type &;
        using const_reference = const Type &;
        using difference_type = cu_uinteger;
        using size_type = cu_uinteger;

    private:
        /** \brief Device data pointer.*/
        value_type *dDataPtr_;

        /** \brief  The size of array elements. */
        size_type size_;

        /** \brief If the pointer is maintained by this. */
        bool maintained_;

        bool zeroCopy_;

    public:
        /**
         * \brief Construct null.
         */
        cu_host inline cu1DArray()
            : dDataPtr_(nullptr), size_(0), maintained_(true), zeroCopy_(false) {}

        /**
         * \brief Construct from size.
         */
        cu_host inline cu1DArray(const size_type size)
            : dDataPtr_(nullptr), size_(size), maintained_(true), zeroCopy_(false) {
            malloc();
        }

        /**
         * \brief Construct from size.
         */
        cu_host inline cu1DArray(const size_type size, const cudaStream_t &streams)
            : dDataPtr_(nullptr), size_(size), maintained_(true), zeroCopy_(false) {
            mallocAsync(streams);
        }

        /**
         * \brief Construct from size and host data pointer.
         */
        cu_host inline cu1DArray(const size_type size, const value_type *__restrict__ hostDataPtr)
            : dDataPtr_(nullptr), size_(size), maintained_(true), zeroCopy_(false) {
            malloc();
            copyFromHost(hostDataPtr);
        }

        /**
         * \brief Construct from size and host data pointer.
         */
        cu_host inline cu1DArray(const size_type size, const value_type *__restrict__ hostDataPtr,
                                 const cudaStream_t &streams)
            : dDataPtr_(nullptr), size_(size), maintained_(true), zeroCopy_(false) {
            mallocAsync(streams);
            copyFromHostAsync(hostDataPtr, streams);
        }

        /**
         * \brief Construct from size and host data pointer and Zero-Copy flag.
         */
        cu_host inline cu1DArray(const size_type size, const value_type *__restrict__ hostDataPtr,
                                 const bool zeroCopy)
            : dDataPtr_(nullptr), size_(size), maintained_(true), zeroCopy_(zeroCopy) {
            if (zeroCopy_) {
                checkCUDAError(
                    cudaHostGetDevicePointer((void **)&dDataPtr_, (void *)hostDataPtr, 0));
            } else {
                malloc();
                copyFromHost(hostDataPtr);
            }
        }

        /**
         * \brief Construct from size and host data pointer and Zero-Copy flag.
         */
        cu_host inline cu1DArray(const size_type size, const value_type *__restrict__ hostDataPtr,
                                 const cudaStream_t &streams, const bool zeroCopy)
            : dDataPtr_(nullptr), size_(size), maintained_(true), zeroCopy_(zeroCopy) {
            if (zeroCopy_) {
                checkCUDAError(
                    cudaHostGetDevicePointer((void **)&dDataPtr_, (void *)hostDataPtr, 0));
            } else {
                mallocAsync(streams);
                copyFromHostAsync(hostDataPtr, streams);
            }
        }

        /**
         * \brief Copy constructor.
         */
        cu_dual inline cu1DArray(const cu1DArray &cc)
            : dDataPtr_(cc.dDataPtr_), size_(cc.size_), maintained_(false),
              zeroCopy_(cc.zeroCopy_) {}

        /**
         * \brief Copy constructor.
         */
        cu_dual inline cu1DArray(cu1DArray &&cc) = delete;

        /**
         * \brief Disallow assignment.
         */
        cu_dual cu1DArray &operator=(const cu1DArray &) = delete;

        /**
         * \brief Disallow assignment.
         */
        cu_dual void operator=(cu1DArray &&cc) = delete;

        /**
         * \brief Destructor.
         */
        cu_dual inline ~cu1DArray() noexcept {}

        /** \brief Device data pointer.*/
        cu_device inline reference operator()(const size_type i) noexcept { return dDataPtr_[i]; }

        /** \brief Device data pointer.*/
        cu_device inline const_reference operator()(const size_type i) const noexcept {
            return dDataPtr_[i];
        }

        /** \brief Device data pointer.*/
        cu_device inline reference operator[](const size_type i) noexcept { return dDataPtr_[i]; }

        /** \brief Device data pointer.*/
        cu_device inline const_reference operator[](const size_type i) const noexcept {
            return dDataPtr_[i];
        }

        /** \brief Device data pointer.*/
        cu_dual inline value_type *dDataPtr() noexcept { return dDataPtr_; }

        /** \brief Device data pointer.*/
        cu_dual inline const value_type *dDataPtr() const noexcept { return dDataPtr_; }

        /** \brief Resize. */
        cu_host inline void resize(const size_type s) {
            if (maintained_) {
                size_ = s;
                if (!zeroCopy_) {
                    malloc();
                }
            }
        }

        /** \brief Clear the array. */
        cu_host inline void clear() {
            if (maintained_) {
                if (!zeroCopy_) {
                    if (dDataPtr_ != nullptr) {
                        checkCUDAError(cudaFree(dDataPtr_));
                        dDataPtr_ = nullptr;
                        size_ = 0;
                    }
                } else {
                    dDataPtr_ = nullptr;
                    size_ = 0;
                }
            }
        }

        /**
         * \brief Return the byte size.
         */
        cu_dual inline size_t byteSize() const noexcept { return size_ * sizeof(value_type); }

        cu_dual inline size_type size() const noexcept { return size_; }

        /** \brief If the pointer is maintained by this. */
        cu_dual inline bool maintained() const noexcept { return maintained_; }

        cu_host inline void copyFromHost(const value_type *__restrict__ hostDataPtr) {
            if (size_ != 0) {
                checkCUDAError(
                    cudaMemcpy(dDataPtr_, hostDataPtr, byteSize(), cudaMemcpyHostToDevice));
            }
        }
        cu_host inline void copyFromHostAsync(const value_type *__restrict__ hostDataPtr,
                                              cudaStream_t stream) {
            if (size_ != 0) {
                checkCUDAError(cudaMemcpyAsync(dDataPtr_, hostDataPtr, byteSize(),
                                               cudaMemcpyHostToDevice, stream));
            }
        }
        cu_host inline void copyToHost(value_type *__restrict__ hostDataPtr) const {
            if (size_ != 0) {
                checkCUDAError(
                    cudaMemcpy(hostDataPtr, dDataPtr_, byteSize(), cudaMemcpyDeviceToHost));
            }
        }
        cu_host inline void copyToHostAsync(value_type *__restrict__ hostDataPtr,
                                            cudaStream_t stream) const {
            if (size_ != 0) {
                checkCUDAError(cudaMemcpyAsync(hostDataPtr, dDataPtr_, byteSize(),
                                               cudaMemcpyDeviceToHost, stream));
            }
        }

    private:
        /**
         * \brief Malloc memory in the device.
         */
        cu_host inline void malloc() {
            if (maintained_) {
                if (size_ != 0) {
                    checkCUDAError(cudaMalloc((void **)&dDataPtr_, byteSize()));
                }
            }
        }

        /**
         * \brief Malloc memory in the device.
         */
        cu_host inline void mallocAsync(const cudaStream_t &streams) {
            if (maintained_) {
                if (size_) {
                    if (cuAttributes::deviceSupportsMemoryPools()) {
                        checkCUDAError(cudaMallocAsync((void **)&dDataPtr_, byteSize(), streams));
                    } else {
                        checkCUDAError(cudaMalloc((void **)&dDataPtr_, byteSize()));
                    }
                }
            }
        }
    };

    /**
     * \brief CUDA 2D array struct.
     */
    template <typename Type> class cu2DArray {
    public:
        using value_type = Type;
        using reference = Type &;
        using const_reference = const Type &;
        using difference_type = cu_uinteger;
        using size_type = cu_uinteger;

    private:
        /** \brief Device data pointer.*/
        value_type *dDataPtr_;

        /** \brief  The returned pitch (or stride) must be used to access array elements. */
        size_t pitch_;

        /** \brief width x height 2D array. */
        size_type width_;

        /** \brief width x height 2D array. */
        size_type height_;

        /** \brief If the pointer is maintained by this. */
        bool maintained_;

    public:
        /**
         * \brief Construct null.
         */
        cu_host inline cu2DArray()
            : dDataPtr_(nullptr), pitch_(size_t(0)), width_(0), height_(0), maintained_(true) {}

        /**
         * \brief Construct from components.
         */
        cu_host inline cu2DArray(const size_type w, const size_type h)
            : dDataPtr_(nullptr), pitch_(size_t(0)), width_(w), height_(h), maintained_(true) {
            malloc();
        }

        /**
         * \brief Construct from components.
         */
        cu_host inline cu2DArray(const size_type w, const size_type h,
                                 const value_type *__restrict__ hostDataPtr)
            : dDataPtr_(nullptr), pitch_(size_t(0)), width_(w), height_(h), maintained_(true) {
            malloc();
            copyFromHost(hostDataPtr);
        }

        /**
         * \brief Disallow copy constructor.
         */
        cu_dual inline cu2DArray(const cu2DArray &cc)
            : dDataPtr_(cc.dDataPtr_), pitch_(cc.pitch_), width_(cc.width_), height_(cc.height_),
              maintained_(false) {}

        /**
         * \brief Copy constructor.
         */
        cu_dual inline cu2DArray(cu2DArray &&cc) = delete;

        /**
         * \brief Disallow assignment.
         */
        cu_dual cu2DArray &operator=(const cu2DArray &) = delete;

        /**
         * \brief Disallow assignment.
         */
        cu_dual void operator=(cu2DArray &&cc) = delete;

        /**
         * \brief Destructor.
         */
        cu_dual inline ~cu2DArray() noexcept {}

        /** \brief Device data pointer.*/
        cu_device inline reference operator()(const size_type i, const size_type j) {
            return ((value_type *)((char *)dDataPtr_ + i * pitch_))[j];
        }

        /** \brief Device data pointer.*/
        cu_device inline const_reference operator()(const size_type i, const size_type j) const {
            return ((value_type *)((char *)dDataPtr_ + i * pitch_))[j];
        }

        cu_device inline value_type *row(const size_type j) {
            return (value_type *)((char *)dDataPtr_ + j * pitch_);
        }

        cu_device inline const value_type *row(const size_type j) const {
            return (value_type *)((char *)dDataPtr_ + j * pitch_);
        }

        /** \brief Device data pointer.*/
        cu_dual inline value_type *dDataPtr() noexcept { return dDataPtr_; }

        /** \brief Device data pointer.*/
        cu_dual inline const value_type *dDataPtr() const noexcept { return dDataPtr_; }

        /** \brief  The returned pitch (or stride) must be used to access array elements. */
        cu_dual inline size_t pitch() const noexcept { return pitch_; }

        /** \brief width x height 2D array. */
        cu_dual inline size_type width() const noexcept { return width_; }

        /** \brief width x height 2D array. */
        cu_dual inline size_type height() const noexcept { return height_; }

        /** \brief Resize. */
        cu_host inline void resize(const size_type w, const size_type h) {
            if (maintained_) {
                width_ = w;
                height_ = h;
                malloc();
            }
        }

        /** \brief Clear the array. */
        cu_host inline void clear() {
            if (maintained_) {
                if (dDataPtr_ != nullptr) {
                    checkCUDAError(cudaFree(dDataPtr_));
                    dDataPtr_ = nullptr;
                    pitch_ = (size_t)0;
                    width_ = 0;
                    height_ = 0;
                }
            }
        }

        /**
         * \brief Return the byty size.
         */
        cu_dual inline size_t byteSize() const noexcept { return width_ * sizeof(value_type); }

        cu_host inline void copyFromHost(const value_type *__restrict__ hostDataPtr) {
            if (width_ != 0 && height_ != 0) {
                checkCUDAError(cudaMemcpy2D(dDataPtr_, pitch_, hostDataPtr, byteSize(), byteSize(),
                                            height_, cudaMemcpyHostToDevice));
            }
        }
        cu_host inline void copyFromHostAsync(const value_type *__restrict__ hostDataPtr,
                                              cudaStream_t stream) {
            if (width_ != 0 && height_ != 0) {
                checkCUDAError(cudaMemcpy2DAsync(dDataPtr_, pitch_, hostDataPtr, byteSize(),
                                                 byteSize(), height_, cudaMemcpyHostToDevice,
                                                 stream));
            }
        }
        cu_host inline void copyToHost(value_type *__restrict__ hostDataPtr) const {
            if (width_ != 0 && height_ != 0) {
                checkCUDAError(cudaMemcpy2D(hostDataPtr, byteSize(), dDataPtr_, pitch_, byteSize(),
                                            height_, cudaMemcpyDeviceToHost));
            }
        }
        cu_host inline void copyToHostAsync(value_type *__restrict__ hostDataPtr,
                                            cudaStream_t stream) const {
            if (width_ != 0 && height_ != 0) {
                checkCUDAError(cudaMemcpy2DAsync(hostDataPtr, byteSize(), dDataPtr_, pitch_,
                                                 byteSize(), height_, cudaMemcpyDeviceToHost,
                                                 stream));
            }
        }

    private:
        /**
         * \brief Malloc memory in the device.
         */
        cu_host inline void malloc() {
            if (maintained_) {
                if (width_ && height_) {
                    checkCUDAError(
                        cudaMallocPitch((void **)&dDataPtr_, &pitch_, byteSize(), height_));
                }
            }
        }

        /**
         * \brief Malloc memory in the device.
         */
        cu_host inline void mallocAsync(cudaStream_t stream) {
            if (maintained_) {
                if (width_ && height_) {
                    checkCUDAError(
                        cudaMallocPitch((void **)&dDataPtr_, &pitch_, byteSize(), height_));
                }
            }
        }
    };

    /**
     * \brief CUDA 3D array struct.
     */
    template <typename Type> class cu3DArray {
    private:
        cudaPitchedPtr dData_;

    public:
        /**
         * \brief Construct null.
         */
        cu_host inline cu3DArray() : dData_() {}

        /**
         * \brief Destructor.
         */
        inline cu_host ~cu3DArray() noexcept {}

        inline cu_host cudaPitchedPtr &dData() noexcept { return dData_; }
        inline cu_host const cudaPitchedPtr &dData() const noexcept { return dData_; }

        /** \brief Clear the array. */
        cu_host inline void clear() { checkCUDAError(cudaFree(dData_.ptr)); }

        /**
         * \brief Make cudaExtent based on input parameters.
         * \param[in] w - Width in elements when referring to array memory
         * \param[in] h - Height in elements
         * \param[in] d - Depth in elements
         */
        inline cu_host cudaExtent makeExtent(const cu_uinteger w, const cu_uinteger h,
                                             const cu_uinteger d) {
            cudaExtent extentCoeff;
            extentCoeff = make_cudaExtent(sizeof(Type) * w, (size_t)h, (size_t)d);
            checkCUDAError(cudaMalloc3D(&dData_, extentCoeff));
            return extentCoeff;
        }

        /**
         * \brief Copy data from host to device.
         * \param[in] hostDataPtr - Data pointer in host memory
         * \param[in] w - Width in elements when referring to array memory
         * \param[in] h - Height in elements
         * \param[in] d - Depth in elements
         */
        cu_host inline void copyFromHost(Type *__restrict__ hostDataPtr, const cu_uinteger w,
                                         const cu_uinteger h, const cu_uinteger d) {
            cudaMemcpy3DParms cpyParmCoeff;
            cpyParmCoeff = {0};
            cpyParmCoeff.srcPtr = make_cudaPitchedPtr((void *)hostDataPtr, sizeof(Type) * w, w, h);
            cpyParmCoeff.dstPtr = dData_;
            cpyParmCoeff.extent = makeExtent(w, h, d);
            cpyParmCoeff.kind = cudaMemcpyHostToDevice;
            checkCUDAError(cudaMemcpy3D(&cpyParmCoeff));
        }

        /**
         * \brief Copy data from host to device.
         * \param[in] hostDataPtr - Data pointer in host memory
         * \param[in] w - Width in elements when referring to array memory
         * \param[in] h - Height in elements
         * \param[in] d - Depth in elements
         * \param[in] stream - Stream identifier
         */
        cu_host inline void copyFromHostAsync(Type *__restrict__ hostDataPtr, const cu_uinteger w,
                                              const cu_uinteger h, const cu_uinteger d,
                                              cudaStream_t stream) {
            cudaMemcpy3DParms cpyParmCoeff;
            cpyParmCoeff = {0};
            cpyParmCoeff.srcPtr = make_cudaPitchedPtr((void *)hostDataPtr, sizeof(Type) * w, w, h);
            cpyParmCoeff.dstPtr = dData_;
            cpyParmCoeff.extent = makeExtent(w, h, d);
            cpyParmCoeff.kind = cudaMemcpyHostToDevice;
            checkCUDAError(cudaMemcpy3DAsync(&cpyParmCoeff, stream));
        }

        /**
         * \brief Copy data from device to host.
         * \param[in] hostDataPtr - Data pointer in host memory
         * \param[in] w - Width in elements when referring to array memory
         * \param[in] h - Height in elements
         * \param[in] d - Depth in elements
         */
        cu_host inline void copyToHost(Type *__restrict__ hostDataPtr, const cu_uinteger w,
                                       const cu_uinteger h, const cu_uinteger d) const {
            cudaMemcpy3DParms cpyParmCoeff;
            cpyParmCoeff = {0};
            cpyParmCoeff.srcPtr = make_cudaPitchedPtr((void *)hostDataPtr, sizeof(Type) * w, w, h);
            cpyParmCoeff.dstPtr = dData_;
            cpyParmCoeff.extent = makeExtent(w, h, d);
            cpyParmCoeff.kind = cudaMemcpyDeviceToHost;
            checkCUDAError(cudaMemcpy3D(&cpyParmCoeff));
        }

        /**
         * \brief Copy data from device to host.
         * \param[in] hostDataPtr - Data pointer in host memory
         * \param[in] w - Width in elements when referring to array memory
         * \param[in] h - Height in elements
         * \param[in] d - Depth in elements
         * \param[in] stream - Stream identifier
         */
        cu_host inline void copyToHostAsync(Type *__restrict__ hostDataPtr, const cu_uinteger w,
                                            const cu_uinteger h, const cu_uinteger d,
                                            cudaStream_t stream) const {
            cudaMemcpy3DParms cpyParmCoeff;
            cpyParmCoeff = {0};
            cpyParmCoeff.srcPtr = make_cudaPitchedPtr((void *)hostDataPtr, sizeof(Type) * w, w, h);
            cpyParmCoeff.dstPtr = dData_;
            cpyParmCoeff.extent = makeExtent(w, h, d);
            cpyParmCoeff.kind = cudaMemcpyDeviceToHost;
            checkCUDAError(cudaMemcpy3DAsync(&cpyParmCoeff, stream));
        }
    };

    class realCu3DArray {
    private:
        cudaPitchedPtr &dData_;

    public:
        /**
         * \brief Construct null.
         */
        cu_dual realCu3DArray() = delete;

        inline cu_dual realCu3DArray(cudaPitchedPtr &dData) : dData_(dData) {}

        cu_dual realCu3DArray(const realCu3DArray &) = delete;

        cu_dual void operator=(const realCu3DArray &) = delete;

        /**
         * \brief Destructor.
         */
        inline cu_dual ~realCu3DArray() noexcept {}

        inline cu_device cu_real &operator()(const cu_uinteger i, const cu_uinteger j,
                                             const cu_uinteger k) {
            char *devPtr = (char *)(dData_.ptr);
            size_t pitch = dData_.pitch;
            size_t slicePitch = pitch * dData_.ysize;
            char *slice = devPtr + i * slicePitch;
            cu_real *row = (cu_real *)(slice + j * pitch);
            return row[k];
        }

        inline cu_device const cu_real &operator()(const cu_uinteger i, const cu_uinteger j,
                                                   const cu_uinteger k) const {
            char *devPtr = (char *)(dData_.ptr);
            size_t pitch = dData_.pitch;
            size_t slicePitch = pitch * dData_.ysize;
            char *slice = devPtr + i * slicePitch;
            cu_real *row = (cu_real *)(slice + j * pitch);
            return row[k];
        }
    };

    class realCu3DArrayConst {
    private:
        const cudaPitchedPtr &dData_;

    public:
        /**
         * \brief Construct null.
         */
        cu_dual realCu3DArrayConst() = delete;

        inline cu_dual realCu3DArrayConst(const cudaPitchedPtr &dData) : dData_(dData) {}

        cu_dual realCu3DArrayConst(const realCu3DArrayConst &) = delete;

        cu_dual void operator=(const realCu3DArrayConst &) = delete;

        /**
         * \brief Destructor.
         */
        inline cu_dual ~realCu3DArrayConst() noexcept {}

        inline cu_device const cu_real &operator()(const cu_uinteger i, const cu_uinteger j,
                                                   const cu_uinteger k) const {
            char *devPtr = (char *)(dData_.ptr);
            size_t pitch = dData_.pitch;
            size_t slicePitch = pitch * dData_.ysize;
            char *slice = devPtr + i * slicePitch;
            cu_real *row = (cu_real *)(slice + j * pitch);
            return row[k];
        }
    };

    /**
     * \brief CUDA managed dynamic array.
     */
    template <class Type, typename sizeType> class cuMArray {
    public:
        using value_type = Type;
        using reference = Type &;
        using const_reference = const Type &;
        using difference_type = sizeType;
        using size_type = sizeType;

    private:
        size_type size_;

        value_type *element_;

        /** \brief If the pointer is maintained by this. */
        bool maintained_;

    public:
        /**
         * \brief Null constructor.
         */
        cu_dual inline cuMArray() : size_(0), element_(nullptr), maintained_(true) {}

        cu_host inline cuMArray(const size_type size)
            : size_(size), element_(nullptr), maintained_(true) {
            malloc();
        }

        cu_host inline cuMArray(const size_type &s, const value_type &elem)
            : size_(s), element_(nullptr), maintained_(true) {
            malloc();
            for (int i = 0; i < size_; ++i) {
                element_[i] = elem;
            }
        }

        /**
         * \brief Copy constructor.
         */
        cu_dual cuMArray(const cuMArray &other)
            : size_(other.size_), element_(other.element_), maintained_(false) {}

        cuMArray &operator=(const cuMArray &other) = delete;

        cu_dual inline ~cuMArray() noexcept {}

        cu_host inline void clear() noexcept {
            size_ = 0;
            free();
        }

    protected:
        /**
         * \brief Free dynamic memory allocation and reserve the size.
         * It is different from clear().
         */
        cu_host inline void free() noexcept {
            if (element_ != nullptr) {
                if (maintained_) {
                    checkCUDAError(cudaFree(element_));
                }
                element_ = nullptr;
            }
        }

        cu_host inline void malloc() {
            if (size_ > 0) {
                if (element_ == nullptr) {
                    checkCUDAError(cudaMallocManaged(&element_, byteSize()));
                }
            }
        }

    public:
        cu_host inline void resize(const size_type newSize) {
            clear();
            maintained_ = true;
            size_ = newSize;
            malloc();
        }

        hur_nodiscard cu_dual inline size_t byteSize() const noexcept {
            return size_ * sizeof(value_type);
        }

        hur_nodiscard cu_dual inline size_type size() const noexcept { return size_; }
        hur_nodiscard cu_dual inline bool empty() const noexcept { return size_ <= 0; }

    public:
        hur_nodiscard cu_dual inline const value_type *data() const noexcept { return element_; }
        hur_nodiscard cu_dual inline value_type *data() noexcept { return element_; }
        hur_nodiscard cu_dual inline const_reference element(const size_type i) const noexcept {
            return element_[i];
        }

        hur_nodiscard cu_dual inline reference element(const size_type i) noexcept {
            return element_[i];
        }

        cu_dual inline void element(value_type &ele, const size_type i) const noexcept {
            ele = element_[i];
        }

        cu_dual inline void replace(const size_type i, const value_type &ele) noexcept {
            element_[i] = ele;
        }
        hur_nodiscard cu_dual inline reference operator[](const size_type i) noexcept {
            return element_[i];
        }

        hur_nodiscard cu_dual inline const_reference operator[](const size_type i) const noexcept {
            return element_[i];
        }

        cu_host inline void memPrefetchAsync(const int devId, cudaStream_t s = 0) const {
            if (size_ > 0) {
                if (element_ != nullptr) {
                    checkCUDAError(cudaMemPrefetchAsync(element_, byteSize(), devId, s));
                }
            }
        }
    };

    /**
     * \brief CUDA managed dynamic block array.
     */
    template <class Type, typename sizeType> class cuMBlockArray {
    public:
        using value_type = Type;

        using size_type = sizeType;

    private:
        size_type size_;

        value_type *element_;

        cu_ushort blockDimx_;
        cu_ushort blockDimy_;

        /** \brief If the pointer is maintained by this. */
        bool maintained_;

    public:
        /**
         * \brief Null constructor.
         */
        cu_dual inline cuMBlockArray() : size_(0), element_(nullptr), maintained_(true) {}

        cu_host inline cuMBlockArray(const size_type size)
            : size_(size), element_(nullptr), blockDimx_(1), blockDimy_(1), maintained_(true) {
            if (!empty()) {
                checkCUDAError(cudaMallocManaged(&element_, byteSize()));
            }
        }

        cu_host inline cuMBlockArray(const size_type size, const cu_ushort blockDim)
            : size_(size), element_(nullptr), blockDimx_(blockDim), blockDimy_(blockDim),
              maintained_(true) {
            if (!empty()) {
                checkCUDAError(cudaMallocManaged(&element_, byteSize()));
            }
        }

        cu_host inline cuMBlockArray(const size_type size, const cu_ushort blockDimx,
                                     const cu_ushort blockDimy)
            : size_(size), element_(nullptr), blockDimx_(blockDimx), blockDimy_(blockDimy),
              maintained_(true) {
            if (!empty()) {
                checkCUDAError(cudaMallocManaged(&element_, byteSize()));
            }
        }

        /**
         * \brief Copy constructor.
         */
        cu_dual cuMBlockArray(const cuMBlockArray &other)
            : size_(other.size_), element_(other.element_), blockDimx_(other.blockDimx_),
              blockDimy_(other.blockDimy_), maintained_(false) {}
        cuMBlockArray &operator=(const cuMBlockArray &other) = delete;

        cu_dual inline ~cuMBlockArray() noexcept {}

        cu_host inline void clear() noexcept {
            size_ = 0;
            blockDimx_ = 0;
            blockDimy_ = 0;
            if (element_ != nullptr) {
                if (maintained_) {
                    checkCUDAError(cudaFree(element_));
                }
                element_ = nullptr;
            }
        }

        cu_host inline void resize(const size_type size, const cu_ushort blockDim) {
            clear();
            size_ = size;
            blockDimx_ = blockDim;
            blockDimy_ = blockDim;
            maintained_ = true;
            if (!empty()) {
                checkCUDAError(cudaMallocManaged(&element_, byteSize()));
            }
        }

        cu_host inline void resize(const size_type size, const cu_ushort blockDimx,
                                   const cu_ushort blockDimy) {
            clear();
            size_ = size;
            blockDimx_ = blockDimx;
            blockDimy_ = blockDimy;
            if (!empty()) {
                checkCUDAError(cudaMallocManaged(&element_, byteSize()));
            }
        }

        hur_nodiscard cu_dual inline size_t byteSize() const noexcept {
            return size_ * blockDimx_ * blockDimy_ * sizeof(value_type);
        }
        hur_nodiscard cu_dual inline cu_ushort blockDimx() const noexcept { return blockDimx_; }
        hur_nodiscard cu_dual inline cu_ushort blockDimy() const noexcept { return blockDimy_; }
        hur_nodiscard cu_dual inline cu_ushort blockSize() const noexcept {
            return blockDimx_ * blockDimy_;
        }
        hur_nodiscard cu_dual inline size_type size() const noexcept { return size_; }
        hur_nodiscard cu_dual inline bool empty() const noexcept {
            return size_ <= 0 || blockDimx_ == 0 || blockDimy_ == 0;
        }

        cu_host inline void memPrefetchAsync(const int devId, cudaStream_t s = 0) const {
            if (size_ > 0) {
                if (element_ != nullptr) {
                    checkCUDAError(cudaMemPrefetchAsync(element_, byteSize(), devId, s));
                }
            }
        }

    public:
        hur_nodiscard cu_dual inline const value_type *data() const noexcept { return element_; }
        hur_nodiscard cu_dual inline value_type *data() noexcept { return element_; }
        class blocks;
        friend class blocks;

        class blocks {
        private:
            cuMBlockArray &arr_;

            size_type id_;

        public:
            blocks() = delete;

            cu_dual inline blocks(cuMBlockArray &arr, const size_type id) : arr_(arr), id_(id) {}

            cu_dual inline blocks(blocks &other) : arr_(other.arr_), id_(other.id_) {}
            cu_dual inline blocks &operator=(const blocks &other) {
                id_ = other.id_;
                return *this;
            }

            hur_nodiscard cu_dual inline value_type *data() noexcept {
                return arr_.data() + id_ * arr_.blockSize();
            }

            hur_nodiscard cu_dual inline const value_type *data() const noexcept {
                return arr_.data() + id_ * arr_.blockSize();
            }

            hur_nodiscard cu_dual inline value_type &operator()(const size_type i,
                                                                const size_type j) noexcept {
                return *(this->data() + i * arr_.blockDimx() + j);
            }

            hur_nodiscard cu_dual inline const value_type &
            operator()(const size_type i, const size_type j) const noexcept {
                return *(this->data() + i * arr_.blockDimx() + j);
            }
        };

        hur_nodiscard cu_dual inline blocks operator()(const size_type id) {
            return blocks(*this, id);
        }

        class constBlocks;
        friend class constBlocks;

        class constBlocks {
        private:
            const cuMBlockArray &arr_;

            size_type id_;

        public:
            constBlocks() = delete;

            cu_dual inline constBlocks(const cuMBlockArray &arr, const size_type id)
                : arr_(arr), id_(id) {}

            cu_dual inline constBlocks(const constBlocks &other)
                : arr_(other.arr_), id_(other.id_) {}
            cu_dual inline constBlocks &operator=(const constBlocks &other) {
                id_ = other.id_;
                return *this;
            }
            cu_dual inline constBlocks &operator=(const blocks &other) {
                id_ = other.id_;
                return *this;
            }

            hur_nodiscard cu_dual inline const value_type *data() const noexcept {
                return arr_.data() + id_ * arr_.blockSize();
            }

            hur_nodiscard cu_dual inline const value_type &
            operator()(const size_type i, const size_type j) const noexcept {
                return *(this->data() + i * arr_.blockDimx() + j);
            }
        };

        hur_nodiscard cu_dual inline constBlocks operator()(const size_type id) const {
            return constBlocks(*this, id);
        }
    };

    /**
     * \brief CUDA dynamic array.
     */
    template <class Type, typename sizeType> class cuArray {
    public:
        using value_type = Type;
        using reference = Type &;
        using const_reference = const Type &;
        using difference_type = sizeType;
        using size_type = sizeType;

    private:
        size_type size_;

        value_type *hElePtr_;
        value_type *dElePtr_;

        /** \brief If the pointer is maintained by this. */
        bool maintained_;

    public:
        /**
         * \brief Null constructor.
         */
        cu_dual inline cuArray()
            : size_(0), hElePtr_(nullptr), dElePtr_(nullptr), maintained_(true) {}

        cu_host inline cuArray(const size_type size)
            : size_(size), hElePtr_(nullptr), dElePtr_(nullptr), maintained_(true) {
            malloc();
        }

        cu_host inline cuArray(const size_type size, const cudaStream_t &s)
            : size_(size), hElePtr_(nullptr), dElePtr_(nullptr), maintained_(true) {
            mallocAsync(s);
        }

        cu_host inline cuArray(const size_type &s, const value_type &elem)
            : size_(s), hElePtr_(nullptr), dElePtr_(nullptr), maintained_(true) {
            malloc();
            for (int i = 0; i < size_; ++i) {
                hElePtr_[i] = elem;
            }
            copyFromHost();
        }

        cu_host inline cuArray(const size_type &s, const value_type *elem)
            : size_(s), hElePtr_(nullptr), dElePtr_(nullptr), maintained_(true) {
            malloc();
            for (int i = 0; i < size_; ++i) {
                hElePtr_[i] = *(elem + i);
            }
            copyFromHost();
        }

        cu_host inline cuArray(const size_type &s, const value_type &elem,
                               const cudaStream_t &stream)
            : size_(s), hElePtr_(nullptr), dElePtr_(nullptr), maintained_(true) {
            mallocAsync(stream);
            for (int i = 0; i < size_; ++i) {
                hElePtr_[i] = elem;
            }
            checkCUDAError(cudaStreamSynchronize(stream));
            copyFromHostAsync(stream);
        }

        cu_host inline cuArray(const size_type &s, const value_type *elem,
                               const cudaStream_t &stream)
            : size_(s), hElePtr_(nullptr), dElePtr_(nullptr), maintained_(true) {
            mallocAsync(stream);
            for (int i = 0; i < size_; ++i) {
                hElePtr_[i] = *(elem + i);
            }
            checkCUDAError(cudaStreamSynchronize(stream));
            copyFromHostAsync(stream);
        }

        /**
         * \brief Copy constructor.
         */
        cu_dual cuArray(const cuArray &other)
            : size_(other.size_), hElePtr_(other.hElePtr_), dElePtr_(other.dElePtr_),
              maintained_(false) {}

        cuArray &operator=(const cuArray &other) = delete;

        cu_dual inline ~cuArray() noexcept {}

        cu_host inline void clear() noexcept {
            size_ = 0;
            free();
        }

    protected:
        /**
         * \brief Free dynamic memory allocation and reserve the size.
         * It is different from clear().
         */
        cu_host inline void free() noexcept {
            if (hElePtr_ != nullptr) {
                if (maintained_) {
                    checkCUDAError(cudaFreeHost(hElePtr_));
                }
                hElePtr_ = nullptr;
            }
            if (dElePtr_ != nullptr) {
                if (maintained_) {
                    checkCUDAError(cudaFree(dElePtr_));
                }
                dElePtr_ = nullptr;
            }
        }

        cu_host inline void malloc() {
            if (maintained_) {
                if (size_ > 0) {
                    if (hElePtr_ == nullptr) {
                        checkCUDAError(
                            cudaHostAlloc((void **)&hElePtr_, byteSize(), cudaHostAllocDefault));
                    }
                    if (dElePtr_ == nullptr) {
                        checkCUDAError(cudaMalloc((void **)&dElePtr_, byteSize()));
                    }
                }
            }
        }

        /**
         * \brief Malloc memory in the device.
         */
        cu_host inline void mallocAsync(const cudaStream_t &streams) {
            if (maintained_) {
                if (size_ > 0) {
                    if (hElePtr_ == nullptr) {
                        checkCUDAError(
                            cudaHostAlloc((void **)&hElePtr_, byteSize(), cudaHostAllocDefault));
                    }

                    if (cuAttributes::deviceSupportsMemoryPools()) {
                        checkCUDAError(cudaMallocAsync((void **)&dElePtr_, byteSize(), streams));
                    } else {
                        checkCUDAError(cudaMalloc((void **)&dElePtr_, byteSize()));
                    }
                }
            }
        }

    public:
        cu_host inline void resize(const size_type newSize) {
            clear();
            maintained_ = true;
            size_ = newSize;
            malloc();
        }

        cu_host inline void resize(const size_type newSize, const cudaStream_t &s) {
            clear();
            maintained_ = true;
            size_ = newSize;
            mallocAsync(s);
        }

        hur_nodiscard cu_dual inline size_t byteSize() const noexcept {
            return size_ * sizeof(value_type);
        }

        hur_nodiscard cu_dual inline size_type size() const noexcept { return size_; }
        hur_nodiscard cu_dual inline bool empty() const noexcept { return size_ <= 0; }

    public:
        hur_nodiscard cu_host inline const value_type *hdata() const noexcept { return hElePtr_; }
        hur_nodiscard cu_host inline value_type *hdata() noexcept { return hElePtr_; }
        hur_nodiscard cu_device inline const value_type *ddata() const noexcept { return dElePtr_; }
        hur_nodiscard cu_device inline value_type *ddata() noexcept { return dElePtr_; }
        hur_nodiscard cu_host inline const_reference helement(const size_type i) const noexcept {
            return hElePtr_[i];
        }

        hur_nodiscard cu_host inline reference helement(const size_type i) noexcept {
            return hElePtr_[i];
        }

        hur_nodiscard cu_device inline const_reference delement(const size_type i) const noexcept {
            return dElePtr_[i];
        }

        hur_nodiscard cu_device inline reference delement(const size_type i) noexcept {
            return dElePtr_[i];
        }

        hur_nodiscard cu_host inline reference operator[](const size_type i) noexcept {
            return hElePtr_[i];
        }

        hur_nodiscard cu_host inline const_reference operator[](const size_type i) const noexcept {
            return hElePtr_[i];
        }

        hur_nodiscard cu_device inline reference operator()(const size_type i) noexcept {
            return dElePtr_[i];
        }

        hur_nodiscard cu_device inline const_reference
        operator()(const size_type i) const noexcept {
            return dElePtr_[i];
        }

        cu_host inline void copyFromHost() {
            if (size_ != 0) {
                checkCUDAError(cudaMemcpy(dElePtr_, hElePtr_, byteSize(), cudaMemcpyHostToDevice));
            }
        }
        cu_host inline void copyFromHostAsync(cudaStream_t stream) {
            if (size_ != 0) {
                checkCUDAError(cudaMemcpyAsync(dElePtr_, hElePtr_, byteSize(),
                                               cudaMemcpyHostToDevice, stream));
            }
        }
        cu_host inline void copyToHost() {
            if (size_ != 0) {
                checkCUDAError(cudaMemcpy(hElePtr_, dElePtr_, byteSize(), cudaMemcpyDeviceToHost));
            }
        }
        cu_host inline void copyToHostAsync(cudaStream_t stream) {
            if (size_ != 0) {
                checkCUDAError(cudaMemcpyAsync(hElePtr_, dElePtr_, byteSize(),
                                               cudaMemcpyDeviceToHost, stream));
            }
        }
    };

    /**
     * \brief CUDA dynamic block array.
     */
    template <class Type, typename sizeType> class cuBlockArray1 {
    public:
        using value_type = Type;

        using size_type = sizeType;

    private:
        size_type size_;

        value_type *hElePtr_;
        value_type *dElePtr_;

        /** \brief  The returned pitch (or stride) must be used to access array elements. */
        size_t pitch_;

        cu_ushort blockDimx_;
        cu_ushort blockDimy_;

        /** \brief If the pointer is maintained by this. */
        bool maintained_;

    public:
        /**
         * \brief Null constructor.
         */
        cu_dual inline cuBlockArray1()
            : size_(0), hElePtr_(nullptr), dElePtr_(nullptr), pitch_(0), blockDimx_(0),
              blockDimy_(0), maintained_(true) {}

        cu_host inline cuBlockArray1(const size_type size)
            : size_(size), hElePtr_(nullptr), dElePtr_(nullptr), pitch_(0), blockDimx_(1),
              blockDimy_(1), maintained_(true) {
            malloc();
        }

        cu_host inline cuBlockArray1(const size_type size, const cu_ushort blockDim)
            : size_(size), hElePtr_(nullptr), dElePtr_(nullptr), pitch_(0), blockDimx_(blockDim),
              blockDimy_(blockDim), maintained_(true) {
            malloc();
        }

        cu_host inline cuBlockArray1(const size_type size, const cu_ushort blockDimx,
                                     const cu_ushort blockDimy)
            : size_(size), hElePtr_(nullptr), dElePtr_(nullptr), pitch_(0), blockDimx_(blockDimx),
              blockDimy_(blockDimy), maintained_(true) {
            malloc();
        }

        /**
         * \brief Copy constructor.
         */
        cu_dual cuBlockArray1(const cuBlockArray1 &other)
            : size_(other.size_), hElePtr_(other.hElePtr_), dElePtr_(other.dElePtr_),
              pitch_(other.pitch_), blockDimx_(other.blockDimx_), blockDimy_(other.blockDimy_),
              maintained_(false) {}
        cuBlockArray1 &operator=(const cuBlockArray1 &other) = delete;

        cu_dual inline ~cuBlockArray1() noexcept {}

        cu_host inline void clear() noexcept {
            pitch_ = 0;
            size_ = 0;
            blockDimx_ = 0;
            blockDimy_ = 0;
            if (hElePtr_ != nullptr) {
                if (maintained_) {
                    checkCUDAError(cudaFreeHost(hElePtr_));
                }
                hElePtr_ = nullptr;
            }

            if (dElePtr_ != nullptr) {
                if (maintained_) {
                    checkCUDAError(cudaFree(dElePtr_));
                }
                dElePtr_ = nullptr;
            }
        }

    protected:
        cu_host inline void malloc() {
            if (maintained_) {
                if (size_ > 0) {
                    if (hElePtr_ == nullptr) {
                        checkCUDAError(cudaHostAlloc((void **)&hElePtr_, size_ * byteSize(),
                                                     cudaHostAllocDefault));
                    }
                    if (dElePtr_ == nullptr) {
                        checkCUDAError(
                            cudaMallocPitch((void **)&dElePtr_, &pitch_, byteSize(), height()));
                    }
                }
            }
        }

    public:
        cu_host inline void resize(const size_type size, const cu_ushort blockDim) {
            clear();
            size_ = size;
            blockDimx_ = blockDim;
            blockDimy_ = blockDim;
            maintained_ = true;
            malloc();
        }

        cu_host inline void resize(const size_type size, const cu_ushort blockDimx,
                                   const cu_ushort blockDimy) {
            clear();
            size_ = size;
            blockDimx_ = blockDimx;
            blockDimy_ = blockDimy;
            malloc();
        }

        hur_nodiscard cu_dual inline size_t byteSize() const noexcept {
            return blockDimx_ * blockDimy_ * sizeof(value_type);
        }
        hur_nodiscard cu_dual inline cu_ushort blockDimx() const noexcept { return blockDimx_; }
        hur_nodiscard cu_dual inline cu_ushort blockDimy() const noexcept { return blockDimy_; }
        hur_nodiscard cu_dual inline cu_ushort blockSize() const noexcept {
            return blockDimx_ * blockDimy_;
        }

        hur_nodiscard cu_dual inline size_type width() const noexcept {
            return blockDimx_ * blockDimy_;
        }

        hur_nodiscard cu_dual inline size_type height() const noexcept { return size_; }

        hur_nodiscard cu_dual inline size_type size() const noexcept { return size_; }
        hur_nodiscard cu_dual inline bool empty() const noexcept {
            return size_ <= 0 || blockDimx_ == 0 || blockDimy_ == 0;
        }

        cu_host inline void copyFromHost() {
            if (!empty()) {
                checkCUDAError(cudaMemcpy2D(dElePtr_, pitch_, hElePtr_, byteSize(), byteSize(),
                                            height(), cudaMemcpyHostToDevice));
            }
        }
        cu_host inline void copyFromHostAsync(cudaStream_t stream) {
            if (!empty()) {
                checkCUDAError(cudaMemcpy2DAsync(dElePtr_, pitch_, hElePtr_, byteSize(), byteSize(),
                                                 height(), cudaMemcpyHostToDevice, stream));
            }
        }
        cu_host inline void copyToHost() {
            if (!empty()) {
                checkCUDAError(cudaMemcpy2D(hElePtr_, byteSize(), dElePtr_, pitch_, byteSize(),
                                            height(), cudaMemcpyDeviceToHost));
            }
        }
        cu_host inline void copyToHostAsync(cudaStream_t stream) {
            if (!empty()) {
                checkCUDAError(cudaMemcpy2DAsync(hElePtr_, byteSize(), dElePtr_, pitch_, byteSize(),
                                                 height(), cudaMemcpyDeviceToHost, stream));
            }
        }

    public:
        hur_nodiscard cu_host inline const value_type *hdata() const noexcept { return hElePtr_; }
        hur_nodiscard cu_host inline value_type *hdata() noexcept { return hElePtr_; }
        hur_nodiscard cu_device inline const value_type *ddata() const noexcept { return dElePtr_; }
        hur_nodiscard cu_device inline value_type *ddata() noexcept { return dElePtr_; }

        class hblocks;
        friend class hblocks;

        class hblocks {
        private:
            cuBlockArray1 &arr_;

            size_type id_;

        public:
            hblocks() = delete;

            cu_host inline hblocks(cuBlockArray1 &arr, const size_type id) : arr_(arr), id_(id) {}

            cu_host inline hblocks(hblocks &other) : arr_(other.arr_), id_(other.id_) {}
            cu_host inline hblocks &operator=(const hblocks &other) {
                id_ = other.id_;
                return *this;
            }

            hur_nodiscard cu_host inline value_type *data() noexcept {
                return arr_.hdata() + id_ * arr_.blockSize();
            }

            hur_nodiscard cu_host inline const value_type *data() const noexcept {
                return arr_.hdata() + id_ * arr_.blockSize();
            }

            hur_nodiscard cu_host inline value_type &operator()(const size_type i,
                                                                const size_type j) noexcept {
                return *(this->data() + i * arr_.blockDimx() + j);
            }

            hur_nodiscard cu_host inline const value_type &
            operator()(const size_type i, const size_type j) const noexcept {
                return *(this->data() + i * arr_.blockDimx() + j);
            }
        };

        hur_nodiscard cu_host inline hblocks operator[](const size_type id) {
            return hblocks(*this, id);
        }

        class dblocks;
        friend class dblocks;

        class dblocks {
        private:
            cuBlockArray1 &arr_;

            size_type id_;

        public:
            dblocks() = delete;

            cu_device inline dblocks(cuBlockArray1 &arr, const size_type id) : arr_(arr), id_(id) {}

            cu_device inline dblocks(dblocks &other) : arr_(other.arr_), id_(other.id_) {}
            cu_device inline dblocks &operator=(const dblocks &other) {
                id_ = other.id_;
                return *this;
            }

            hur_nodiscard cu_device inline value_type *data() noexcept {
                return (value_type *)((char *)arr_.ddata() + id_ * arr_.pitch_);
            }

            hur_nodiscard cu_device inline const value_type *data() const noexcept {
                return (value_type *)((char *)arr_.ddata() + id_ * arr_.pitch_);
            }

            hur_nodiscard cu_device inline value_type &operator()(const size_type i,
                                                                  const size_type j) noexcept {
                return *(this->data() + i * arr_.blockDimx() + j);
            }

            hur_nodiscard cu_device inline const value_type &
            operator()(const size_type i, const size_type j) const noexcept {
                return *(this->data() + i * arr_.blockDimx() + j);
            }
        };

        hur_nodiscard cu_device inline dblocks operator()(const size_type id) {
            return dblocks(*this, id);
        }

        class consthBlocks;
        friend class consthBlocks;

        class consthBlocks {
        private:
            const cuBlockArray1 &arr_;

            size_type id_;

        public:
            consthBlocks() = delete;

            cu_host inline consthBlocks(const cuBlockArray1 &arr, const size_type id)
                : arr_(arr), id_(id) {}

            cu_host inline consthBlocks(const consthBlocks &other)
                : arr_(other.arr_), id_(other.id_) {}
            cu_host inline consthBlocks &operator=(const consthBlocks &other) {
                id_ = other.id_;
                return *this;
            }
            cu_host inline consthBlocks &operator=(const hblocks &other) {
                id_ = other.id_;
                return *this;
            }

            hur_nodiscard cu_host inline const value_type *data() const noexcept {
                return arr_.data() + id_ * arr_.blockSize();
            }

            hur_nodiscard cu_host inline const value_type &
            operator()(const size_type i, const size_type j) const noexcept {
                return *(this->data() + i * arr_.blockDimx() + j);
            }
        };

        hur_nodiscard cu_host inline consthBlocks operator[](const size_type id) const {
            return consthBlocks(*this, id);
        }

        class constdblocks;
        friend class constdblocks;

        class constdblocks {
        private:
            const cuBlockArray1 &arr_;

            size_type id_;

        public:
            constdblocks() = delete;

            cu_device inline constdblocks(const cuBlockArray1 &arr, const size_type id)
                : arr_(arr), id_(id) {}

            cu_device inline constdblocks(constdblocks &other) : arr_(other.arr_), id_(other.id_) {}
            cu_device inline constdblocks &operator=(const constdblocks &other) {
                id_ = other.id_;
                return *this;
            }

            cu_device inline constdblocks &operator=(const dblocks &other) {
                id_ = other.id_;
                return *this;
            }

            hur_nodiscard cu_device inline const value_type *data() const noexcept {
                return (value_type *)((char *)arr_.ddata() + id_ * arr_.pitch_);
            }

            hur_nodiscard cu_device inline const value_type &
            operator()(const size_type i, const size_type j) const noexcept {
                return *(this->data() + i * arr_.blockDimx() + j);
            }
        };

        hur_nodiscard cu_device inline constdblocks operator()(const size_type id) const {
            return constdblocks(*this, id);
        }
    };

    /**
     * \brief CUDA dynamic block array.
     */
    template <class Type, typename sizeType> class cuBlockArray2 {
    public:
        using value_type = Type;

        using size_type = sizeType;

    private:
        size_type size_;

        value_type *hElePtr_;
        value_type *dElePtr_;

        cu_ushort blockDimx_;
        cu_ushort blockDimy_;

        /** \brief If the pointer is maintained by this. */
        bool maintained_;

    public:
        /**
         * \brief Null constructor.
         */
        cu_dual inline cuBlockArray2()
            : size_(0), hElePtr_(nullptr), dElePtr_(nullptr), blockDimx_(0), blockDimy_(0),
              maintained_(true) {}

        cu_host inline cuBlockArray2(const size_type size)
            : size_(size), hElePtr_(nullptr), dElePtr_(nullptr), blockDimx_(1), blockDimy_(1),
              maintained_(true) {
            malloc();
        }

        cu_host inline cuBlockArray2(const size_type size, const cu_ushort blockDim)
            : size_(size), hElePtr_(nullptr), dElePtr_(nullptr), blockDimx_(blockDim),
              blockDimy_(blockDim), maintained_(true) {
            malloc();
        }

        cu_host inline cuBlockArray2(const size_type size, const cu_ushort blockDimx,
                                     const cu_ushort blockDimy)
            : size_(size), hElePtr_(nullptr), dElePtr_(nullptr), blockDimx_(blockDimx),
              blockDimy_(blockDimy), maintained_(true) {
            malloc();
        }

        /**
         * \brief Copy constructor.
         */
        cu_dual cuBlockArray2(const cuBlockArray2 &other)
            : size_(other.size_), hElePtr_(other.hElePtr_), dElePtr_(other.dElePtr_),
              blockDimx_(other.blockDimx_), blockDimy_(other.blockDimy_), maintained_(false) {}
        cuBlockArray2 &operator=(const cuBlockArray2 &other) = delete;

        cu_dual inline ~cuBlockArray2() noexcept {}

        cu_host inline void clear() noexcept {
            size_ = 0;
            blockDimx_ = 0;
            blockDimy_ = 0;
            if (hElePtr_ != nullptr) {
                if (maintained_) {
                    checkCUDAError(cudaFreeHost(hElePtr_));
                }
                hElePtr_ = nullptr;
            }

            if (dElePtr_ != nullptr) {
                if (maintained_) {
                    checkCUDAError(cudaFree(dElePtr_));
                }
                dElePtr_ = nullptr;
            }
        }

    protected:
        cu_host inline void malloc() {
            if (maintained_) {
                if (size_ > 0) {
                    if (hElePtr_ == nullptr) {
                        checkCUDAError(
                            cudaHostAlloc((void **)&hElePtr_, byteSize(), cudaHostAllocDefault));
                    }
                    if (dElePtr_ == nullptr) {
                        checkCUDAError(cudaMalloc((void **)&dElePtr_, byteSize()));
                    }
                }
            }
        }

    public:
        cu_host inline void resize(const size_type size, const cu_ushort blockDim) {
            clear();
            size_ = size;
            blockDimx_ = blockDim;
            blockDimy_ = blockDim;
            maintained_ = true;
            malloc();
        }

        cu_host inline void resize(const size_type size, const cu_ushort blockDimx,
                                   const cu_ushort blockDimy) {
            clear();
            size_ = size;
            blockDimx_ = blockDimx;
            blockDimy_ = blockDimy;
            malloc();
        }

        hur_nodiscard cu_dual inline size_t byteSize() const noexcept {
            return size_t(size_ * blockDimx_ * blockDimy_) * sizeof(value_type);
        }
        hur_nodiscard cu_dual inline cu_ushort blockDimx() const noexcept { return blockDimx_; }
        hur_nodiscard cu_dual inline cu_ushort blockDimy() const noexcept { return blockDimy_; }
        hur_nodiscard cu_dual inline cu_ushort blockSize() const noexcept {
            return blockDimx_ * blockDimy_;
        }

        hur_nodiscard cu_dual inline size_type size() const noexcept { return size_; }

        hur_nodiscard cu_dual inline bool empty() const noexcept {
            return size_ <= 0 || blockDimx_ == 0 || blockDimy_ == 0;
        }

        cu_host inline void copyFromHost() {
            if (!empty()) {
                checkCUDAError(cudaMemcpy(dElePtr_, hElePtr_, byteSize(), cudaMemcpyHostToDevice));
            }
        }
        cu_host inline void copyFromHostAsync(cudaStream_t stream) {
            if (!empty()) {
                checkCUDAError(cudaMemcpyAsync(dElePtr_, hElePtr_, byteSize(),
                                               cudaMemcpyHostToDevice, stream));
            }
        }
        cu_host inline void copyToHost() {
            if (!empty()) {
                checkCUDAError(cudaMemcpy(hElePtr_, dElePtr_, byteSize(), cudaMemcpyDeviceToHost));
            }
        }
        cu_host inline void copyToHostAsync(cudaStream_t stream) {
            if (!empty()) {
                checkCUDAError(cudaMemcpyAsync(hElePtr_, dElePtr_, byteSize(),
                                               cudaMemcpyDeviceToHost, stream));
            }
        }

    public:
        hur_nodiscard cu_host inline const value_type *hdata() const noexcept { return hElePtr_; }
        hur_nodiscard cu_host inline value_type *hdata() noexcept { return hElePtr_; }
        hur_nodiscard cu_device inline const value_type *ddata() const noexcept { return dElePtr_; }
        hur_nodiscard cu_device inline value_type *ddata() noexcept { return dElePtr_; }

        class hblocks;
        friend class hblocks;

        class hblocks {
        private:
            cuBlockArray2 &arr_;

            size_type id_;

        public:
            hblocks() = delete;

            cu_host inline hblocks(cuBlockArray2 &arr, const size_type id) : arr_(arr), id_(id) {}

            cu_host inline hblocks(hblocks &other) : arr_(other.arr_), id_(other.id_) {}
            cu_host inline hblocks &operator=(const hblocks &other) {
                id_ = other.id_;
                return *this;
            }

            hur_nodiscard cu_host inline value_type *data() noexcept {
                return arr_.hdata() + id_ * arr_.blockSize();
            }

            hur_nodiscard cu_host inline const value_type *data() const noexcept {
                return arr_.hdata() + id_ * arr_.blockSize();
            }

            hur_nodiscard cu_host inline value_type &operator()(const size_type i,
                                                                const size_type j) noexcept {
                return *(this->data() + i * arr_.blockDimx() + j);
            }

            hur_nodiscard cu_host inline const value_type &
            operator()(const size_type i, const size_type j) const noexcept {
                return *(this->data() + i * arr_.blockDimx() + j);
            }
        };

        hur_nodiscard cu_host inline hblocks operator[](const size_type id) {
            return hblocks(*this, id);
        }

        class dblocks;
        friend class dblocks;

        class dblocks {
        private:
            cuBlockArray2 &arr_;

            size_type id_;

        public:
            dblocks() = delete;

            cu_device inline dblocks(cuBlockArray2 &arr, const size_type id) : arr_(arr), id_(id) {}

            cu_device inline dblocks(dblocks &other) : arr_(other.arr_), id_(other.id_) {}
            cu_device inline dblocks &operator=(const dblocks &other) {
                id_ = other.id_;
                return *this;
            }

            hur_nodiscard cu_device inline value_type *data() noexcept {
                return arr_.ddata() + id_ * arr_.blockSize();
            }

            hur_nodiscard cu_device inline const value_type *data() const noexcept {
                return arr_.ddata() + id_ * arr_.blockSize();
            }

            hur_nodiscard cu_device inline value_type &operator()(const size_type i,
                                                                  const size_type j) noexcept {
                return *(this->data() + i * arr_.blockDimx() + j);
            }

            hur_nodiscard cu_device inline const value_type &
            operator()(const size_type i, const size_type j) const noexcept {
                return *(this->data() + i * arr_.blockDimx() + j);
            }
        };

        hur_nodiscard cu_device inline dblocks operator()(const size_type id) {
            return dblocks(*this, id);
        }

        class consthBlocks;
        friend class consthBlocks;

        class consthBlocks {
        private:
            const cuBlockArray2 &arr_;

            size_type id_;

        public:
            consthBlocks() = delete;

            cu_host inline consthBlocks(const cuBlockArray2 &arr, const size_type id)
                : arr_(arr), id_(id) {}

            cu_host inline consthBlocks(const consthBlocks &other)
                : arr_(other.arr_), id_(other.id_) {}
            cu_host inline consthBlocks &operator=(const consthBlocks &other) {
                id_ = other.id_;
                return *this;
            }
            cu_host inline consthBlocks &operator=(const hblocks &other) {
                id_ = other.id_;
                return *this;
            }

            hur_nodiscard cu_host inline const value_type *data() const noexcept {
                return arr_.data() + id_ * arr_.blockSize();
            }

            hur_nodiscard cu_host inline const value_type &
            operator()(const size_type i, const size_type j) const noexcept {
                return *(this->data() + i * arr_.blockDimx() + j);
            }
        };

        hur_nodiscard cu_host inline consthBlocks operator[](const size_type id) const {
            return consthBlocks(*this, id);
        }

        class constdblocks;
        friend class constdblocks;

        class constdblocks {
        private:
            const cuBlockArray2 &arr_;

            size_type id_;

        public:
            constdblocks() = delete;

            cu_device inline constdblocks(const cuBlockArray2 &arr, const size_type id)
                : arr_(arr), id_(id) {}

            cu_device inline constdblocks(constdblocks &other) : arr_(other.arr_), id_(other.id_) {}
            cu_device inline constdblocks &operator=(const constdblocks &other) {
                id_ = other.id_;
                return *this;
            }

            cu_device inline constdblocks &operator=(const dblocks &other) {
                id_ = other.id_;
                return *this;
            }

            hur_nodiscard cu_device inline const value_type *data() const noexcept {
                return arr_.ddata() + id_ * arr_.blockSize();
            }

            hur_nodiscard cu_device inline const value_type &
            operator()(const size_type i, const size_type j) const noexcept {
                return *(this->data() + i * arr_.blockDimx() + j);
            }
        };

        hur_nodiscard cu_device inline constdblocks operator()(const size_type id) const {
            return constdblocks(*this, id);
        }
    };

#ifndef destroyCuArray
// Destroy arrays stored in GPU gloabl memory.
#define destroyCuArray(arr) arr.clear()
#endif // !destroyCuArray

} // namespace OpenHurricane
#endif // CUDA_PARALLEL