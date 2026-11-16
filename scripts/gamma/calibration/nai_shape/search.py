import numpy as np
from . import model


def anchors_to_p(da, db, ya, yb, n_b):
    """Преобразует сдвиги в двух опорных точках в вектор параметров p = [d, gN, 0, 0].

    Параметры
    ---------
    da, db : float
        Сдвиги (каналы) в точках ya и yb.
    ya, yb : float
        Координаты опорных точек на шкале фона (ya < yb).
    n_b : float
        Общее число каналов фона.

    Возвращает
    ----------
    numpy.ndarray, shape (4,)
        Вектор параметров [d, gN, 0, 0].
    """
    gN = (db - da) * n_b / (yb - ya)
    d = da - gN * ya / n_b
    return np.array([d, gN, 0.0, 0.0])


def anchor_grid_fit(prob, free, n_b, ya, yb, fwa, fwb,
                    span_a, span_b, step_fwhm=0.25, n_starts=5):
    """Глобальный поиск начального приближения по сетке сдвигов в опорных точках.

    По двум развязанным сдвигам δa и δb (в точках ya, yb шкалы фона)
    строится сетка χ², находятся локальные минимумы, из каждого запускается
    локальный спуск (model.fit), выбирается лучший результат.

    Параметры
    ---------
    prob : object
        Объект с методом `chi2(p)` для вычисления χ².
    free : list of int
        Индексы свободных параметров в p = [d, gN, q, s].
    n_b : float
        Общее число каналов фона.
    ya, yb : float
        Координаты опорных точек (ya < yb).
    fwa, fwb : float
        ПШПВ (FWHM) в каналах в точках ya и yb.
    span_a, span_b : float
        Полуширины поиска сдвига (каналы).
    step_fwhm : float, default=0.25
        Шаг сетки в долях FWHM.
    n_starts : int, default=5
        Число лучших кандидатов, из которых запускается локальный спуск.

    Возвращает
    ----------
    dict
        Ключи: "p", "chi2", "ndf", "minima", "grid", "da_grid", "db_grid".
    """
    # Построение сеток сдвигов
    da_grid = np.arange(-span_a, span_a + 1e-9, step_fwhm * fwa)
    db_grid = np.arange(-span_b, span_b + 1e-9, step_fwhm * fwb)

    na = len(da_grid)
    nb = len(db_grid)

    # Заполнение матрицы χ²
    G = np.empty((na, nb), dtype=np.float64)
    for i in range(na):
        for j in range(nb):
            p = anchors_to_p(da_grid[i], db_grid[j], ya, yb, n_b)
            G[i, j] = prob.chi2(p)

    # Поиск локальных минимумов в окрестности 3×3 (точка ≤ всех соседей)
    candidates = []
    for i in range(na):
        for j in range(nb):
            is_min = True
            for di in range(-1, 2):
                for dj in range(-1, 2):
                    if di == 0 and dj == 0:
                        continue
                    ni, nj = i + di, j + dj
                    if 0 <= ni < na and 0 <= nj < nb:
                        if G[i, j] > G[ni, nj]:
                            is_min = False
                            break
                if not is_min:
                    break
            if is_min:
                candidates.append((G[i, j], da_grid[i], db_grid[j]))

    # Сортировка по возрастанию χ², отбор n_starts лучших
    candidates.sort(key=lambda x: x[0])
    candidates = candidates[:n_starts]

    # Локальный спуск из каждого кандидата (без сеток)
    results = []
    for chi2_start, da, db in candidates:
        p0 = anchors_to_p(da, db, ya, yb, n_b)
        r = model.fit(prob, free, p0)
        results.append({
            "p": r["p"].tolist(),
            "chi2": float(r["chi2"]),
            "start": [float(da), float(db)],
        })

    # Сортировка результатов по возрастанию χ²
    results.sort(key=lambda x: x["chi2"])

    # Лучший результат
    best = results[0] if results else None
    if best is None:
        # Сеточный минимум не найден — используем минимум всей сетки как запасной
        i_min, j_min = np.unravel_index(np.argmin(G), G.shape)
        p0 = anchors_to_p(da_grid[i_min], db_grid[j_min], ya, yb, n_b)
        r = model.fit(prob, free, p0)
        results.append({
            "p": r["p"].tolist(),
            "chi2": float(r["chi2"]),
            "start": [float(da_grid[i_min]), float(db_grid[j_min])],   # ручная правка: было da_grid[j_min] (VERIFY.md п. 7)
        })
        results.sort(key=lambda x: x["chi2"])
        best = results[0]
        p_best = np.array(best["p"], dtype=np.float64)
        chi2_best = best["chi2"]
        ndf_best = r["ndf"]
    else:
        # Перезапускаем model.fit для лучшего, чтобы получить ndf
        i_best, j_best = None, None
        # Находим, какой кандидат дал лучший результат (по start)
        for idx, res in enumerate(results):
            if idx == 0:
                da_v, db_v = res["start"]
                break
        p_best = anchors_to_p(da_v, db_v, ya, yb, n_b)
        r_best = model.fit(prob, free, p_best)
        chi2_best = float(r_best["chi2"])
        ndf_best = r_best["ndf"]
        p_best = r_best["p"]

    return {
        "p": p_best,
        "chi2": chi2_best,
        "ndf": ndf_best,
        "minima": results,
        "grid": G,
        "da_grid": da_grid.tolist(),
        "db_grid": db_grid.tolist(),
    }
