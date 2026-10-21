"""World map data for the spread map: country outlines from Natural Earth 110m (public domain, `resources/world_110m.json`) and a resolver from the country names found in
genome metadata (NCBI / Pathoplexus spellings, old names, small territories the 110m map leaves out) to a map entry."""
import json
import unicodedata
from pathlib import Path

_FILE = Path(__file__).parent / "resources" / "world_110m.json"
_CACHE: dict = {}

# spellings that differ from Natural Earth's names
ALIASES = {
    "usa": "united states of america", "united states": "united states of america", "us": "united states of america", "uk": "united kingdom", "great britain": "united kingdom",
    "viet nam": "vietnam", "cote d'ivoire": "ivory coast", "cote divoire": "ivory coast", "zaire": "democratic republic of the congo", "drc": "democratic republic of the congo",
    "dr congo": "democratic republic of the congo", "congo, democratic republic of the": "democratic republic of the congo", "democratic republic of congo": "democratic republic of the congo",
    "republic of the congo": "republic of the congo", "congo": "republic of the congo", "russian federation": "russia", "korea, republic of": "south korea", "republic of korea": "south korea",
    "korea, south": "south korea", "korea, democratic people's republic of": "north korea", "czech republic": "czechia", "swaziland": "eswatini", "turkey": "turkey", "turkiye": "turkey",
    "tanzania, united republic of": "united republic of tanzania", "tanzania": "united republic of tanzania", "laos": "laos", "lao pdr": "laos", "syria": "syria", "iran, islamic republic of": "iran",
    "bolivia, plurinational state of": "bolivia", "venezuela, bolivarian republic of": "venezuela", "taiwan, province of china": "taiwan", "macedonia": "north macedonia", "burma": "myanmar",
    "the gambia": "gambia", "timor-leste": "east timor", "cabo verde": "cape verde", "bahamas, the": "the bahamas", "bahamas": "the bahamas", "serbia": "republic of serbia",
}

# places the 110m outlines do not have: (longitude, latitude)
POINTS = {
    "singapore": ("SGP", "Singapore", 103.82, 1.35), "hong kong": ("HKG", "Hong Kong", 114.17, 22.32), "macao": ("MAC", "Macao", 113.54, 22.2), "malta": ("MLT", "Malta", 14.38, 35.94),
    "bahrain": ("BHR", "Bahrain", 50.55, 26.07), "mauritius": ("MUS", "Mauritius", 57.55, -20.25), "reunion": ("REU", "Réunion", 55.54, -21.12), "seychelles": ("SYC", "Seychelles", 55.49, -4.68),
    "maldives": ("MDV", "Maldives", 73.5, 3.2), "barbados": ("BRB", "Barbados", -59.54, 13.19), "martinique": ("MTQ", "Martinique", -61.02, 14.64), "guadeloupe": ("GLP", "Guadeloupe", -61.55, 16.25),
    "puerto rico": ("PRI", "Puerto Rico", -66.5, 18.22), "cayman islands": ("CYM", "Cayman Islands", -81.25, 19.3), "bermuda": ("BMU", "Bermuda", -64.75, 32.3), "monaco": ("MCO", "Monaco", 7.42, 43.74),
    "luxembourg": ("LUX", "Luxembourg", 6.13, 49.82), "andorra": ("AND", "Andorra", 1.6, 42.5), "liechtenstein": ("LIE", "Liechtenstein", 9.55, 47.17), "gibraltar": ("GIB", "Gibraltar", -5.35, 36.14),
    "french polynesia": ("PYF", "French Polynesia", -149.4, -17.68), "new caledonia": ("NCL", "New Caledonia", 165.6, -21.25), "mayotte": ("MYT", "Mayotte", 45.17, -12.83),
    "sao tome and principe": ("STP", "São Tomé and Príncipe", 6.61, 0.19), "comoros": ("COM", "Comoros", 43.87, -11.88), "cape verde": ("CPV", "Cabo Verde", -23.6, 15.1),
    "faroe islands": ("FRO", "Faroe Islands", -6.9, 61.9), "curacao": ("CUW", "Curaçao", -68.99, 12.17), "aruba": ("ABW", "Aruba", -69.97, 12.52), "guam": ("GUM", "Guam", 144.79, 13.44),
}


POINT_CONTINENT = {"ABW": "North America", "AND": "Europe", "BHR": "Asia", "BMU": "North America", "BRB": "North America", "COM": "Africa", "CPV": "Africa", "CUW": "North America", "CYM": "North America",
                   "FRO": "Europe", "GIB": "Europe", "GLP": "North America", "GUM": "Oceania", "HKG": "Asia", "LIE": "Europe", "LUX": "Europe", "MAC": "Asia", "MCO": "Europe", "MDV": "Asia",
                   "MLT": "Europe", "MTQ": "North America", "MUS": "Africa", "MYT": "Africa", "NCL": "Oceania", "PRI": "North America", "PYF": "Oceania", "REU": "Africa", "SGP": "Asia",
                   "STP": "Africa", "SYC": "Africa"}


def norm(s: str) -> str:
    s = unicodedata.normalize("NFKD", (s or "").strip()).encode("ascii", "ignore").decode().lower()
    return " ".join(s.split())


def world() -> dict:
    if not _CACHE:
        data = json.loads(_FILE.read_text(encoding="utf-8"))
        idx = {}
        for c in data["countries"]:
            for k in (c["name"], c["long"], c["admin"]):
                idx.setdefault(norm(k), c)
        _CACHE["data"], _CACHE["idx"] = data, idx
    return _CACHE["data"]


def resolve(name: str) -> dict | None:
    """{iso, name, lon, lat, continent} of a country name as written in genome metadata, or None (blank, 'missing', unknown)."""
    world()
    n = norm(name)
    if not n or n in ("missing", "unknown", "none", "na", "n/a", "not applicable", "not collected", "not provided"):
        return None
    n = ALIASES.get(n, n)
    c = _CACHE["idx"].get(n)
    if c:
        return {"iso": c["iso"], "name": c["name"], "lon": c["c"][0], "lat": c["c"][1], "continent": c.get("cont", "")}
    p = POINTS.get(n)
    if p:
        return {"iso": p[0], "name": p[1], "lon": p[2], "lat": p[3], "continent": POINT_CONTINENT.get(p[0], "")}
    return None


def outlines() -> bytes:
    """The outlines file as served to the browser."""
    return _FILE.read_bytes()
