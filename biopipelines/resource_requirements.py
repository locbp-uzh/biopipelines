# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""What a tool needs from the machine, and the warning when a batch's request does not meet it.

A tool declares ``RESOURCE_REQUIREMENTS`` as a class attribute. Only the ``gpu`` key is read::

    RESOURCE_REQUIREMENTS = {"gpu": {"incompatible": ["H100"], "use": "A100"}}

``incompatible`` lists the GPU models the tool does not run on. ``compatible``, when not empty, lists the only models
it runs on, so any model outside it counts as incompatible. ``use`` is the request the warning suggests.

A batch's ``Resources(gpu=...)`` is judged against that declaration when the batch scripts are generated:

* a request that names a model is judged by its name alone;
* a request that does not pin a model (``any``, ``gpu``, a count, ``high-memory``, a memory such as ``80GB``,
  ``!X``, ``A|B``) is judged against ``machine.gpu_models`` of the config, the GPU models the machine offers, and
  warns when it may land on an incompatible one. A model whose memory the config leaves empty can be the answer to
  any memory request. Without ``machine.gpu_models`` only a request that names a model is judged;
* ``none``, an empty request and a token that is neither a model nor a memory are not judged.

The result is a list of messages, never an error: the scheduler decides where the job lands, and the tool may run.
"""

from typing import Dict, List, Optional, Sequence, Tuple

GPU_KEYS = ("compatible", "incompatible", "use")

# Models the Resources() documentation names; a request for one of them is a model, not an unknown token.
KNOWN_GPU_MODELS = ("T4", "L4", "V100", "A100", "H100", "H200")

# The memory tokens the SLURM backend turns into a GPUMEM constraint, and what "high-memory" stands for.
MEMORY_TOKENS = ("24GB", "32GB", "80GB", "96GB", "140GB")
HIGH_MEMORY_TOKENS = ("32GB", "80GB", "96GB")

ANY_REQUESTS = ("any", "gpu")
NOT_JUDGED = ("", "none")


def validate_requirements(tool_name: str, requirements) -> None:
    """Refuse a declaration the check would silently ignore: a misspelled key, a value that is not a list of names."""
    if not isinstance(requirements, dict):
        raise ValueError(f"{tool_name}: RESOURCE_REQUIREMENTS must be a dict, got {type(requirements).__name__}.")
    gpu = requirements.get("gpu")
    if gpu is None:
        return
    if not isinstance(gpu, dict):
        raise ValueError(f"{tool_name}: RESOURCE_REQUIREMENTS['gpu'] must be a dict with {GPU_KEYS}.")
    unknown = sorted(set(gpu) - set(GPU_KEYS))
    if unknown:
        raise ValueError(f"{tool_name}: RESOURCE_REQUIREMENTS['gpu'] has no key {unknown}; the keys are {GPU_KEYS}.")
    for key in ("compatible", "incompatible"):
        names = gpu.get(key, [])
        if not isinstance(names, (list, tuple)) or not all(isinstance(name, str) and name for name in names):
            raise ValueError(f"{tool_name}: RESOURCE_REQUIREMENTS['gpu']['{key}'] must be a list of GPU model names.")
    if "use" in gpu and not (isinstance(gpu["use"], str) and gpu["use"]):
        raise ValueError(f"{tool_name}: RESOURCE_REQUIREMENTS['gpu']['use'] must be a GPU request such as 'A100'.")


def _models_with_memory(tokens: Sequence[str], inventory: Dict[str, str]) -> List[str]:
    """Inventory models a memory request can be answered by; a model with no memory stated can be any of them."""
    wanted = {token.upper() for token in tokens}
    return [model for model, memory in inventory.items() if not memory or str(memory).upper() in wanted]


def _is_model(token: str, inventory: Dict[str, str]) -> bool:
    upper = token.upper()
    return upper in KNOWN_GPU_MODELS or upper in {model.upper() for model in inventory}


def resolve_request(request, inventory: Dict[str, str]) -> Tuple[List[str], List[str]]:
    """(named, possible): the models a GPU request names, and the inventory models it may land on without naming one."""
    if request is None:
        return [], []
    text = str(request).strip()
    lower = text.lower()
    if lower in NOT_JUDGED:
        return [], []
    if lower in ANY_REQUESTS or lower.isdigit():
        return [], list(inventory)
    if lower == "high-memory":
        return [], _models_with_memory(HIGH_MEMORY_TOKENS, inventory)
    if text.startswith("!"):
        excluded = text[1:].strip().upper()
        return [], [model for model in inventory if model.upper() != excluded]
    named: List[str] = []
    possible: List[str] = []
    for token in (part.strip() for part in text.split("|")):
        if token.upper() in MEMORY_TOKENS:
            possible.extend(_models_with_memory([token], inventory))
        elif _is_model(token, inventory):
            named.append(token)
    return named, possible


def gpu_warnings(tool_name: str, request, requirements: Optional[dict], inventory: Optional[Dict[str, str]] = None) -> List[str]:
    """The messages for one tool's declaration against one batch's GPU request; empty when it can run there."""
    gpu = (requirements or {}).get("gpu") or {}
    incompatible = {name.upper() for name in gpu.get("incompatible", [])}
    compatible = {name.upper() for name in gpu.get("compatible", [])}
    if not incompatible and not compatible:
        return []

    def unsupported(model: str) -> bool:
        return model.upper() in incompatible or (bool(compatible) and model.upper() not in compatible)

    named, possible = resolve_request(request, inventory or {})
    fix = gpu.get("use") or (sorted(gpu.get("compatible", []))[0] if compatible else None)
    advice = f' Request Resources(gpu="{fix}") instead.' if fix else ""

    messages: List[str] = []
    refused = [model for model in dict.fromkeys(named) if unsupported(model)]
    if refused:
        messages.append(
            f'{tool_name} does not support {", ".join(refused)}, which Resources(gpu="{request}") requests.{advice}')
    risky = [model for model in dict.fromkeys(possible) if unsupported(model) and model not in refused]
    if risky:
        messages.append(
            f'{tool_name}: Resources(gpu="{request}") may land on {", ".join(risky)}, '
            f'which {tool_name} does not support.{advice}')
    return messages
