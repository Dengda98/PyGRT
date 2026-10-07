// 在 RTD 的 Settings → Addons → Custom Script 中配置此文件的在线地址
// 使用 /zh-cn/latest/_static/vercount_rtd.js，确保旧版本也加载最新脚本
(function () {
    "use strict";

    function loadVercount() {
        // 最新版本已在 footer.html 中加载 Vercount，避免重复计数
        if (document.querySelector('script[src="https://events.vercount.one/js"]')) {
            return;
        }

        // 将旧版本的显示标签交给 Vercount，避免不蒜子继续更新这些标签
        ["site_pv", "site_uv", "page_pv"].forEach(function (metric) {
            ["value", "container"].forEach(function (kind) {
                var oldId = "busuanzi_" + kind + "_" + metric;
                var newId = "vercount_" + kind + "_" + metric;
                var element = document.getElementById(oldId);

                if (element && !document.getElementById(newId)) {
                    element.id = newId;
                }
            });
        });

        var script = document.createElement("script");
        script.src = "https://events.vercount.one/js";
        script.async = true;
        document.head.appendChild(script);
    }

    if (document.readyState === "loading") {
        document.addEventListener("DOMContentLoaded", loadVercount, { once: true });
    } else {
        loadVercount();
    }
})();
