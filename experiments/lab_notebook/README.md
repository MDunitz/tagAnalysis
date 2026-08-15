## `lab_notebook`

Lab notebook entries are organized by year. Each entry follows the naming convention:

```
lab_notebook/
├── YYYYMMDD_PROTOCOL_NAME_TEMPLATE.md   # Template for new entries
├── 2025/
│   └── 20250618_BIOFILM_GROWTH.md       # Example entry
├── 2026/
│   └── ...
```

### Creating a new entry

1. Copy `YYYYMMDD_PROTOCOL_NAME_TEMPLATE.md` into the appropriate year directory (create the year folder if it doesn't exist yet, e.g., `2027/`)
2. Rename using the date the experiment starts: `YYYYMMDD_PROTOCOL_NAME.md`
3. Fill in sampling date, who performed the protocol, and all measurements
4. When complete, update the `status` and `reason` in the frontmatter
5. **Never copy an existing entry and rename it** — always start from the template

### Year folders

Year folders are created on demand. When you add the first entry of a new year, create the directory (e.g., `2027/`) and place the entry inside it — no `.gitkeep` placeholder is needed, because the directory is non-empty as soon as the entry lands.
