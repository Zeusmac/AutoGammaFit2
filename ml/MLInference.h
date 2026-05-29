#pragma once

#ifdef HAS_ONNX

#include <onnxruntime_cxx_api.h>
#include <array>
#include <string>
#include <vector>

// Shared ONNX Runtime environment — one instance per process.
inline Ort::Env& OrtGlobalEnv()
{
    static Ort::Env env(ORT_LOGGING_LEVEL_WARNING, "AutoGammaFit");
    return env;
}

// Thin wrapper around one ONNX session.
// Load once at startup; call RunScalar() once per sliding window.
struct OrtMLModel {
    Ort::Session session_;

    explicit OrtMLModel(const std::string& path)
        : session_(OrtGlobalEnv(), path.c_str(), Ort::SessionOptions{})
    {}

    // Run the model on one window; returns the single scalar output.
    float RunScalar(const std::vector<float>& window)
    {
        Ort::AllocatorWithDefaultOptions alloc;
        auto in_name  = session_.GetInputNameAllocated(0, alloc);
        auto out_name = session_.GetOutputNameAllocated(0, alloc);

        std::array<int64_t, 2> shape{1, static_cast<int64_t>(window.size())};
        auto mem = Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);
        auto in_tensor = Ort::Value::CreateTensor<float>(
            mem,
            const_cast<float*>(window.data()), window.size(),
            shape.data(), shape.size());

        const char* in_names[]  = {in_name.get()};
        const char* out_names[] = {out_name.get()};
        auto results = session_.Run(Ort::RunOptions{nullptr},
                                    in_names,  &in_tensor, 1,
                                    out_names, 1);
        return *results[0].GetTensorData<float>();
    }
};

#endif // HAS_ONNX
