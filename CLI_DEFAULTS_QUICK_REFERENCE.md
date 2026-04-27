# CLI Defaults Quick Reference Card

## 🐛 Bugs Fixed in 3DRegBsplines

| Parameter | Old | New | Issue | Impact |
|-----------|-----|-----|-------|--------|
| `--lambda` | 0 | 0.5 | NGF disabled | ❌→✅ |
| `--nu` | 0 | 0 | MSE off (OK) | 📝 Documented |
| `--costfunctionconvergencefactor` | 1e12 | 1e7 | Too loose | 🔴→✅ |
| `--overlappadding` | 1 | 5 | Too small | ⚠️→✅ |

---

## Quick Test

```bash
cd /data/PROJECTS/registrationSuite
rm -rf build && mkdir build && cd build
cmake .. && make -j$(nproc)

# Verify fix
./bin/3DRegBsplines --help | grep "lambda\|convergence\|padding"
```

Expected: `lambda 0.5`, `convergence 1e7`, `padding 5`

---

## Affected File
- `src/3DRegistration/3DRegBsplines/src/3DRegBsplines.cxx`
  - Lines: 94, 100, 104, 153

---

## Documentation Files
1. `CLI_DEFAULTS_AUDIT.md` — Full analysis
2. `FIXES_APPLIED.md` — Changes & rebuild
3. `CLI_DEFAULTS_VERIFICATION.md` — Complete summary
4. `CLI_DEFAULTS_QUICK_REFERENCE.md` — This file

