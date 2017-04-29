#pragma once

#include "spdlog/spdlog.h"

/* Define global logging macros */
#define INFO(msg, ...) spdlog::get("LOG")->info(msg, ##__VA_ARGS__);
#define WARNING(msg, ...) spdlog::get("LOG")->warn(msg, ##__VA_ARGS__);
#define ERROR(msg, ...) spdlog::get("LOG")->error(msg, ##__VA_ARGS__);
