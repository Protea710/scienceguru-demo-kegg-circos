# KEGG ortholog circos plot — reproducible demo

一张三层结构的环状图，全部数值来自两个公开数据源，可用一条命令重画。

![figure](output/figure.png)

## 30 秒复现

![30 秒复现](assets/repro.gif)

干净克隆后执行 `bash run.sh`：前 21 秒下载两个数据源，随后 8 秒建表、出图、跑完 129 项校验。
上面这段是**真实执行的终端输出按 1.2× 回放**（下载段已跳过），不是演示动画——每一行的数字都来自那次运行。

## 这张图在说什么

| 层 | 视觉元素 | 编码的数据 |
|---|---|---|
| 最外圈热图 | 3 圈色块，每圈一个样本 | 48 个 KEGG ortholog 在 LA1 / LB1 / LC1 三个样本中的相对丰度（%） |
| 内圈连线 | 黄 / 蓝曲线 | 每个 KO 连到它所属的 KEGG 功能类；线宽 = 该 KO 参与的 KEGG 通路数（connection degree，0–16） |
| 右侧弧 | 淡色轨道 + 3 根柱 | 每个功能类被选中的 KO 丰度之和；外圈三条细弧是样本配色 |

黄色 = KEGG 一级类「Metabolism」，蓝色 = 其余四个一级类。

## 数据来源

1. **功能类归属**：KEGG BRITE `ko00001`（KO → pathway class 层级），`https://rest.kegg.jp/get/br:ko00001`。
   脚本每次运行时重新抓取，抓取时间与 sha256 记录在 `data/provenance.json`。
   KEGG 数据的版权归 Kanehisa Laboratories，学术用途。
2. **KO 丰度**：figshare `10.6084/m9.figshare.29645648.v1`，标题 *KEGG KO abundance_Mollisols MR SR CT*，
   授权 **CC BY 4.0**。原始表 94 列 × 9437 个 KO，样本名 = 采样点字母 + 处理字母 + 重复号
   （A = 单施化肥对照，B = 施有机肥，C = 秸秆还田）。

本图使用同一采样点 L 的三个处理：**LA1（control）、LB1（manure）、LC1（straw）**。
没有任何模拟值：`data/provenance.json` 里 `simulated_values` 为 `false`。

## 版式与可读性（第二版重做）

画布按 1.2:1 重新定比例，标题两级、关键文字加粗，配色对比度实测如下（WCAG 相对亮度比，白色底）：

| 元素 | 颜色 | 对比度 |
|---|---|---|
| 功能类轨道 · Metabolism | `#D6E4F1` 填充 / `#6F94B8` 描边 | 1.29 : 1 / 3.18 : 1 |
| 功能类轨道 · 其他 | `#F2CFDC` 填充 / `#B96A82` 描边 | 1.43 : 1 / 3.86 : 1 |
| 热图中点 | `#F1EBED` | 1.18 : 1（原为 `#FAF4F6`，1.09 : 1，几乎不可见） |

第一版的淡粉填充只有 1.17 : 1，在白底上等于没有形状，现已改到 1.43 : 1 并加 3.9 : 1 的描边。

## 复现

```bash
bash run.sh          # 抓数据 → 建表 → 出图 → 校验
```

单独运行：

```bash
python3 code/01_prepare_data.py   # 下载 + 解析 + 生成三张中间表
python3 code/02_plot_circos.py    # 渲染 output/figure.{png,svg,pdf}
python3 code/03_verify.py         # 129 项一致性检查，失败即非零退出
```

依赖：Python 3.9+、matplotlib（`pip install -r requirements.txt`）。**不需要 R、不需要 circlize、不需要 GPU、不联网也能重画**
（第一次运行会下载两个原始文件，之后离线可重复）。

## 数据处理的两条规则

1. **相对丰度**：每个样本内，以「落在图中 17 个功能类里的 4138 个 KO」的计数总和为分母，换算成百分比。
2. **KO 选取**：每个功能类取均值丰度最高的 4 个 KO，合并去重后按丰度截取 48 个，最后按 KO 编号排序
   （对应原图按编号排布的观感）。规则写在 `01_prepare_data.py` 的常量里，改一行即可换一种选法。

## 校验做了什么

`code/03_verify.py` 共 129 项断言，全部通过才算成功，包括：

- 原始文件 sha256 与 `provenance.json` 记录一致（数据没被悄悄替换）
- 48 个 KO 都能在 KEGG BRITE 与 figshare 原始表中同时找到
- `ko_links.tsv` 的每一行、每个功能类的 `n_links`、每个 KO 的 `classes` 列互相吻合
- 每个功能类的三样本合计 = 其成员 KO 丰度之和（误差 < 1e-4）
- 三份图像文件确实写出且非空

## 与原参考图的差异

原图是一份宏基因组结果的截图，没有公开原始矩阵，因此这里**不是像素级临摹**，而是在真实数据上重建同一种图式。有意保留的差异：

- KO 编号改为沿半径方向排列。原图沿弧线切向排列，48 个编号在切向会互相压字，径向可读性更好；
  编号带与热图环之间留了 0.018 个坐标单位的空隙，脚本每次渲染都会打印实测的最外侧半径做自检。
- 功能类标签与 KO 编号统一为细体，只比编号大 1pt（4.8 → 5.8），单层放射排布。早先用 8.2pt 粗体、
  两层交错，是因为粗字在窄扇区里必然压字；换成细字后角间距足够，交错反而成了多余的不规则。
- 样本配色改为深蓝 / 琥珀 / 紫，避开热图的青—粉渐变，否则「粉色柱子」和「粉色高值格子」会被读成同一个意思。
- connection degree 标尺从弧线改为左下角直线刻度，避免与连线重叠。
- 颜色标尺按 5%–95% 分位截断并把中点设在样本中位数，否则 0.13%–1.12% 的长尾会让绝大多数格子接近白色；
  截断事实写在图上，不做隐藏。
- 画布由程序按内容包围盒自动拟合，不靠手工调坐标；改字号或字数后不会出现标签被裁掉或大片空白。
- 画布比例、标题文案、字号、配色都是 `code/02_plot_circos.py` 顶部的常量，换比例只需改 `FIG_W / FIG_H`。
- 底部图例块整体按 `LEG_SCALE` 缩放（尺寸、字号、间距同时乘它），当前取 0.78。
- 功能类扇形的弧长按「该类被选中的 KO 数」分配，柱高按「这些 KO 的丰度之和」，两者都只描述图中呈现的子集，不代表该类在完整数据中的总量。

## 文件

```
code/01_prepare_data.py   抓取、解析、生成 data/processed/*.tsv
code/02_plot_circos.py    单文件渲染器（matplotlib，386 行）
code/03_verify.py         一致性校验
data/raw/                 两个原始文件（KEGG 快照、figshare 表）
data/processed/           ko_selected.tsv / ko_links.tsv / class_summary.tsv
data/provenance.json      来源、授权、sha256、推导说明
output/figure.png|svg|pdf
run.sh                    一键复现

说明：写作用稿与录屏分镜保存在本地 `docs/`，未随本仓库发布。
```

## 发布到 git 时的注意事项

- figshare 表是 CC BY 4.0，可以随仓库分发，请保留出处与作者署名。
- 代码与图采用 MIT 许可（见 `LICENSE`）；上游数据各自保留原授权：figshare 表为 CC BY 4.0，
  KEGG 数据版权归 Kanehisa Laboratories、限学术用途，因此未随仓库分发。
- KEGG 快照建议**不要**提交：让 `run.sh` 重新抓取，既不违反 KEGG 的分发条款，也正好证明脚本能独立取数。
  若 KEGG 更新了 `ko00001`，`03_verify.py` 会因为 sha256 不符而报警——那是数据源变了，不是图错了。
